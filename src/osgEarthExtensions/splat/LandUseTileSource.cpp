/* -*-c++-*- */
/* osgEarth - Dynamic map generation toolkit for OpenSceneGraph
 * Copyright 2008-2014 Pelican Mapping
 * http://osgearth.org
 *
 * osgEarth is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include "LandUseTileSource"
#include <osgEarth/ImageLayer>
#include <osgEarth/MapFrame>
#include <osgEarth/Registry>
#include <osgEarth/ImageUtils>
#include <osgEarthUtil/SimplexNoise>

using namespace osgEarth;
using namespace osgEarth::Splat;

namespace
{
    osg::Vec2 getSplatCoords(const TileKey& key, float baseLOD, const osg::Vec2& covUV)
    {
        osg::Vec2 out;

        float dL = (float)key.getLOD() - baseLOD;
        float factor = pow(2.0f, dL);
        float invFactor = 1.0/factor;
        out.set( covUV.x()*invFactor, covUV.y()*invFactor ); 

        // For upsampling we need to calculate an offset as well
        if ( factor >= 1.0 )
        {
            unsigned wide, high;
            key.getProfile()->getNumTiles(key.getLOD(), wide, high);

            float tileX = (float)key.getTileX();
            float tileY = (float)(wide-1-key.getTileY()); // swap Y. (not done in the shader version.)

            osg::Vec2 a( floor(tileX*invFactor), floor(tileY*invFactor) );
            osg::Vec2 b( a.x()*factor, a.y()*factor );
            osg::Vec2 c( (a.x()+1.0f)*factor, (a.y()+1.0f)*factor );
            osg::Vec2 offset( (tileX-b.x())/(c.x()-b.x()), (tileY-b.y())/(c.y()-b.y()) );

            out += offset;
        }

        return out;
    }

    osg::Vec2 warpCoverageCoords(const osg::Vec2& covIn, float noise, float warp)
    {
        float n1 = 2.0 * noise - 1.0;
        return osg::Vec2(
            osg::clampBetween( covIn.x() + n1*warp, 0.0f, 1.0f ),
            osg::clampBetween( covIn.y() + n1*warp, 0.0f, 1.0f ) );
    }

    float getNoise(osgEarth::Util::SimplexNoise& noiseGen, const osg::Vec2& uv)
    {
        // TODO: check that u and v are 0..s and not 0..s-1
        double n = noiseGen.getTiledValue(uv.x(), uv.y());
        n = osg::clampBetween(n, 0.0, 1.0);
        //out = n;
        return n;
    }
}


LandUseTileSource::LandUseTileSource(const LandUseOptions& options) :
TileSource( options ),
_options  ( options )
{
    //nop
}

TileSource::Status
LandUseTileSource::initialize(const osgDB::Options* dbOptions)
{
    _dbOptions = Registry::instance()->cloneOrCreateOptions(dbOptions);

    const Profile* profile = getProfile();
    if ( !profile )
    {
        profile = osgEarth::Registry::instance()->getGlobalGeodeticProfile();
        setProfile( profile );
    }

    // load the image layer:
    if ( _options.imageLayerOptions().isSet() )
    {
        ImageLayerOptions ilo = _options.imageLayerOptions().get();
        ilo.cachePolicy() = CachePolicy::NO_CACHE;
        _imageLayer = new ImageLayer( ilo );
        _imageLayer->setTargetProfileHint( profile );
    }

    // set up the IO options so that we do not cache input data.
    CachePolicy::NO_CACHE.apply( _dbOptions.get() );

    // set up the noise generator.
    const float F[4] = { 4.0f, 16.0f, 4.0f, 8.0f };
    const float P[4] = { 0.8f,  0.6f, 0.8f, 0.9f };
    const float L[4] = { 2.2f,  1.7f, 3.0f, 4.0f };
    
    // Configure the noise function:
    _noiseGen.setNormalize  ( true );
    _noiseGen.setRange      ( 0.0, 1.0 );
    _noiseGen.setFrequency  ( F[0] );
    _noiseGen.setPersistence( P[0] );
    _noiseGen.setLacunarity ( L[0] );
    _noiseGen.setOctaves    ( 8 );

    return STATUS_OK;
}

osg::Image*
LandUseTileSource::createImage(const TileKey&    key,
                               ProgressCallback* progress)
{
    if ( !_imageLayer.valid() )
        return 0L;

    // fetch the image for this key, using a fallback loop until we get data.
    GeoImage image;
    for(TileKey k = key; k.valid() && !image.valid(); k = k.createParentKey())
    {
        image = _imageLayer->createImage(k, progress);
    }

    if (!image.valid() )
        return 0L;
    
    // calculate the subwindow of our key inside the source image, which will be non-unit
    // if we had to fall back.
    float scale = key.getExtent().width() / image.getExtent().width();
    osg::Vec2 bias;
    bias.x() = (key.getExtent().xMin() - image.getExtent().xMin()) / image.getExtent().width();
    bias.y() = (key.getExtent().yMin() - image.getExtent().yMin()) / image.getExtent().height();

    osg::Image* in = image.getImage();
    osg::Image* out = new osg::Image();

    // Allocate a suitable format:
    GLenum type;
    GLint  internalFormat;

    if ( _options.bits().isSetTo(8u) )
    {
        // 8-bit integer:
        type           = GL_UNSIGNED_BYTE;
        internalFormat = GL_LUMINANCE8;
    }
    else if ( _options.bits().isSetTo(16u) )
    {
        // 16-bit integer:
        type           = GL_UNSIGNED_SHORT;
        internalFormat = GL_LUMINANCE16;
    }
    else
    {
        // 32-bit float:
        type           = GL_FLOAT;
        internalFormat = GL_LUMINANCE32F_ARB;
    }
    
    out->allocateImage(256, 256, 1, GL_LUMINANCE, type);
    out->setInternalTextureFormat(internalFormat);

    float noiseLOD = _options.baseLOD().get();
    float warp     = _options.warpFactor().get();

    osg::Vec2 cov;    // coverage coordinates
    float     noise;  // noise value
    osg::Vec2 noiseCoords;

    float du = 1.0f / (float)(in->s()-1);
    float dv = 1.0f / (float)(in->t()-1);

    ImageUtils::PixelReader read ( in );
    ImageUtils::PixelWriter write( out );

    for(float u=0.0f; u<=1.0f; u+=du)
    {
        for(float v=0.0f; v<=1.0f; v+=dv)
        {
            osg::Vec2 cov(scale*u + bias.x(), scale*v + bias.y());

            // Noise is like a repeating overlay at the noiseLOD. So sample it using
            // straight U/V tile coordinates.
            noiseCoords = getSplatCoords( key, noiseLOD, osg::Vec2(u,v) );
            noise = getNoise( _noiseGen, noiseCoords );

            cov = warpCoverageCoords(cov, noise, warp);

            osg::Vec4 texel = read(cov.x(), cov.y());
            write.f(texel, u, v);

            //testing: visualize noise
            //write.f(osg::Vec4f(noise,noise,noise,1), u, v);
        }
    }

    return out;
}
