/* -*-c++-*- */
/* osgEarth - Geospatial SDK for OpenSceneGraph
 * Copyright 2020 Pelican Mapping
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
#include "WindLayer"
#include <osgEarth/Shaders>
#include <osgEarth/StringUtils>
#include <osgEarth/GeoTransform>
#include <osg/Texture3D>
#include <osg/Program>
#include <osg/BindImageTexture>
#include <osgUtil/CullVisitor>

using namespace osgEarth;

#define LC "[WindLayer] "

REGISTER_OSGEARTH_LAYER(wind, WindLayer);

#define WIND_DIM_X 8
#define WIND_DIM_Y 8
#define WIND_DIM_Z 16
#define RADIUS 75

//........................................................................

namespace 
{
    struct WindData // keep me 16-byte aligned
    {
        GLfloat position[4];
        GLfloat direction[3];
        GLfloat speed;
    };

    struct DrawState
    {
        DrawState() :
            _buffer(INT_MAX),
            _bufferSize(0) { }

        GLuint _buffer;
        osg::GLintptr _bufferSize;
    };

    struct CameraState
    {
        CameraState() :
            _windData(NULL), 
            _numWindsAllocated(0) { }

        WindData* _windData;
        unsigned _numWindsAllocated;

        osg::ref_ptr<osg::Uniform> _viewToTexMatrix;
        osg::ref_ptr<osg::Uniform> _texToViewMatrix;
        osg::ref_ptr<osg::StateSet> _stateSet;
        osg::ref_ptr<osg::Texture3D> _texture;
    };

    struct WindDrawable : public osg::Drawable
    {
    public:
        WindDrawable();

        void setup(const osgDB::Options* readOptions);

        void drawImplementation(osg::RenderInfo& ri) const;

        void compileGLObjects(osg::RenderInfo& ri) const;
        void releaseGLObjects(osg::State* state) const;
        void resizeGLObjectBuffers(unsigned maxSize);
        void updateBuffers(const osg::Camera* camera);

        osg::ref_ptr<const osgDB::Options> _readOptions;

        std::vector<osg::ref_ptr<Wind> > _winds;

        mutable osg::buffered_object<DrawState> _ds;

        mutable CameraState _cs;
    };

    WindDrawable::WindDrawable()
    {
        setCullingActive(false);
        setDataVariance(osg::Object::DYNAMIC); // so we can update the wind instances synchronously
    }

    void WindDrawable::setup(const osgDB::Options* readOptions)
    {
        //TODO: make this per-camera state
        osg::StateSet* computeSS = getOrCreateStateSet();

        // Install the compute shader that will generate the texture
        Shaders shaders;
        std::string source = ShaderLoader::load(shaders.WindComputer, shaders, _readOptions.get());
        osg::Shader* computeShader = new osg::Shader(osg::Shader::COMPUTE, source);
        osg::Program* program = new osg::Program();
        program->addShader(computeShader);
        computeSS->setAttribute(program, 1);

        // Make our wind texture
        // TODO: move this to the per-camera stateset
        CameraState& cs = _cs;

        cs._texture = new osg::Texture3D();
        cs._texture->setTextureSize(WIND_DIM_X, WIND_DIM_Y, WIND_DIM_Z);
        cs._texture->setInternalFormat(GL_RGBA8);
        cs._texture->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR);
        cs._texture->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
        cs._texture->setWrap(osg::Texture::WRAP_S, osg::Texture::CLAMP_TO_EDGE);
        cs._texture->setWrap(osg::Texture::WRAP_T, osg::Texture::CLAMP_TO_EDGE);
        cs._texture->setWrap(osg::Texture::WRAP_R, osg::Texture::CLAMP_TO_EDGE);

        // binding so the compute shader can write to the texture
        computeSS->addUniform(new osg::Uniform("oe_wind_tex", 0));
        computeSS->setAttribute(new osg::BindImageTexture(0, cs._texture.get(), osg::BindImageTexture::WRITE_ONLY, GL_RGBA8, 0, GL_TRUE));

        //_matrixUniform = new osg::Uniform("oe_wind_matrix", osg::Matrixf::identity());
        //ss->addUniform(_matrixUniform.get()); // remove this...
        // no add; done in the sharedstateset

        cs._texToViewMatrix = new osg::Uniform("oe_wind_texToViewMatrix", osg::Matrixf::identity());
        computeSS->addUniform(cs._texToViewMatrix.get());

        computeSS->setRenderBinDetails(-90210, "RenderBin");
    }

    void WindDrawable::compileGLObjects(osg::RenderInfo& ri) const
    {
        osg::State* state = ri.getState();
        if (state)
        {
            osg::GLExtensions* ext = state->get<osg::GLExtensions>();

            DrawState& ds = _ds[state->getContextID()];
            CameraState& cs = _cs;

            GLuint requiredBufferSize = sizeof(WindData) * (_winds.size()+1);

            if (ds._buffer == INT_MAX || ds._bufferSize < requiredBufferSize)
            {
                if (ds._buffer != INT_MAX)
                    ext->glDeleteBuffers(1, &ds._buffer);

                ext->glGenBuffers(1, &ds._buffer);

                ds._bufferSize = requiredBufferSize;

                ext->glBindBuffer(GL_SHADER_STORAGE_BUFFER, ds._buffer);
                ext->glBufferStorage(GL_SHADER_STORAGE_BUFFER, ds._bufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
            }

            // download to GPU
            ext->glBindBuffer(GL_SHADER_STORAGE_BUFFER, ds._buffer);
            ext->glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, ds._bufferSize, cs._windData);
        }
    }

    void WindDrawable::updateBuffers(const osg::Camera* camera)
    {
        CameraState& cs = _cs;

        if (cs._numWindsAllocated < _winds.size()+1)
        {
            if (cs._windData != NULL)
                delete [] cs._windData;

            // add one for the terminator.
            cs._windData = new WindData[_winds.size()+1];
            ::memset(cs._windData, 0, sizeof(WindData)*(_winds.size()+1));

            cs._numWindsAllocated = _winds.size()+1;
        }

        size_t i;
        for(i=0; i<_winds.size(); ++i)
        {
            const Wind* wind = _winds[i].get();

            if (wind->type() == Wind::TYPE_POINT)
            {
                // transform from world to camera-view space
                osg::Vec3d posView =  wind->getPointWorld() * camera->getViewMatrix();

                cs._windData[i].position[0] = posView.x();
                cs._windData[i].position[1] = posView.y();
                cs._windData[i].position[2] = posView.z();
                cs._windData[i].position[3] = 1.0;
            }
            else // TYPE_DIRECTIONAL
            {
                // transform from world to camera-view space
                osg::Vec3f dir;
                dir.x() = wind->direction()->x();
                dir.y() = wind->direction()->y();
                dir = osg::Matrixf::transform3x3(dir, camera->getViewMatrix());
                dir.normalize();

                cs._windData[i].direction[0] = dir.x();
                cs._windData[i].direction[1] = dir.y();
                cs._windData[i].direction[2] = dir.z();
                cs._windData[i].position[3] = 0.0f;
            }

            cs._windData[i].speed = wind->speed()->as(Units::KNOTS);
        }

        // terminator record: set power to negative.
        cs._windData[i].speed = -1.0f;
    }

    void WindDrawable::releaseGLObjects(osg::State* state) const
    {
        //TODO
    }

    void WindDrawable::resizeGLObjectBuffers(unsigned maxSize)
    {
        _ds.resize(maxSize);
    }

    void WindDrawable::drawImplementation(osg::RenderInfo& ri) const
    {
        DrawState& ds = _ds[ri.getState()->getContextID()];
        osg::GLExtensions* ext = ri.getState()->get<osg::GLExtensions>();

        // update buffer with wind data
        compileGLObjects(ri);

        // activate layout() binding point:
        ext->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ds._buffer);

        // run it
        ext->glDispatchCompute(WIND_DIM_X, WIND_DIM_Y, WIND_DIM_Z);

        // sync the output
        ext->glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
}

//........................................................................

// default: gentle west-to-east breeze
Wind::Wind() :
    _type(TYPE_DIRECTIONAL),
    _direction(osg::Vec2f(1,0)),
    _speed(Speed(5.0, Units::KNOTS)),
    _pointWorld(osg::Vec3d(0,0,0))
{
    //nop
}

void
Wind::setPoint(const GeoPoint& point)
{
    _point = point;
    point.toWorld(_pointWorld);
}

Wind::Wind(const Config& conf)
{
    conf.get("type", "point", type(), TYPE_POINT);
    conf.get("type", "directional", type(), TYPE_DIRECTIONAL);
    conf.get("point", point());
    conf.get("direction", direction());
    conf.get("speed", speed());
}

Config
Wind::getConfig() const
{
    Config conf("wind");
    conf.set("type", "point", type(), TYPE_POINT);
    conf.set("type", "directional", type(), TYPE_DIRECTIONAL);
    conf.set("point", point());
    conf.set("direction", direction());
    conf.set("speed", speed());
    return conf;
}

//........................................................................

Config
WindLayer::Options::getConfig() const
{
    Config conf = Layer::Options::getConfig();
    if (!winds().empty())
    {
        Config windsConf("winds");
        for(std::vector<osg::ref_ptr<Wind> >::const_iterator i = winds().begin();
            i != winds().end();
            ++i)
        {
            windsConf.add("wind", i->get()->getConfig());
        }
        conf.add(windsConf);
    }
    return conf;
}

void
WindLayer::Options::fromConfig(const Config& conf)
{
    ortho().setDefault(true);
    winds().clear();

    const ConfigSet windsConf = conf.child("winds").children();
    for(ConfigSet::const_iterator i = windsConf.begin(); i != windsConf.end(); ++i)
    {
        Wind* wind = new Wind(*i);
        winds().push_back(wind);
    }
}

//........................................................................

void
WindLayer::addWind(Wind* wind)
{
    WindDrawable* wd = static_cast<WindDrawable*>(_drawable.get());
    wd->_winds.push_back(wind);
}

void
WindLayer::init()
{
    Layer::init();

    // Never cache decals
    layerHints().cachePolicy() = CachePolicy::NO_CACHE;

    _radius = RADIUS;
}

osg::Node*
WindLayer::getNode() const
{
    return _drawable.get();
}

void
WindLayer::setTerrainResources(TerrainResources* res)
{
    res->reserveTextureImageUnit(_unit, "WindLayer");

    // Create the wind drawable that will provide a wind texture
    WindDrawable* wd = new WindDrawable();
    wd->setup(getReadOptions());
    _drawable = wd;

#if 0 // TESTING
    Wind* wind = new Wind();
    wind->type() = Wind::TYPE_DIRECTIONAL;
    wind->direction()->set(-1.0f, 0.0f);
    wind->speed()->set(10, Units::KNOTS);
    addWind(wind);

    Wind* wind2 = new Wind();
    wind2->type() = Wind::TYPE_DIRECTIONAL;
    wind2->direction()->set(0.0f, -1.0f);
    wind2->speed()->set(10, Units::KNOTS);
    addWind(wind2);
#endif

    for(int i=0; i<options().winds().size(); ++i)
    {
        addWind(options().winds()[i].get());
    }
}

osg::StateSet*
WindLayer::getSharedStateSet(osg::NodeVisitor* nv) const
{
    osgUtil::CullVisitor* cv = static_cast<osgUtil::CullVisitor*>(nv);

    WindDrawable* windDrawable = static_cast<WindDrawable*>(_drawable.get());

    CameraState& cs = windDrawable->_cs;

    //todo: per camera
    if (!cs._stateSet.valid())
    {
        cs._stateSet = new osg::StateSet();

        cs._viewToTexMatrix = new osg::Uniform("oe_wind_matrix", osg::Matrixf::identity());
        cs._stateSet->addUniform(cs._viewToTexMatrix.get());
        cs._stateSet->setDefine("OE_WIND_TEX_MATRIX", "oe_wind_matrix");

        cs._stateSet->addUniform(new osg::Uniform("oe_wind_tex", _unit.unit()));
        cs._stateSet->setTextureAttribute(_unit.unit(), cs._texture.get(), osg::StateAttribute::ON);
        cs._stateSet->setDefine("OE_WIND_TEX", "oe_wind_tex");
    }

    // this xforms from clip [-1..1] to texture [0..1] space
    static osg::Matrix clipToTexture = 
        osg::Matrix::translate(1.0,1.0,1.0) * 
        osg::Matrix::scale(0.5,0.5,0.5);

    osg::Matrix rttProjection;

    if (options().ortho() == true)
    {
        rttProjection = osg::Matrix::ortho(
            -_radius, _radius,
            -_radius, _radius,
            0, _radius);
    }
    else
    { 
        double y,a,n,f;
        cv->getCurrentCamera()->getProjectionMatrix().getPerspective(y,a,n,f);

        // pushing the NEAR out reduces jitter a lot
        rttProjection.makePerspective(y, a, 5.0, _radius);
    }

    // view to texture:
    osg::Matrix camViewToTexture = rttProjection * clipToTexture;
    cs._viewToTexMatrix->set(camViewToTexture);

    // texture to view (for the compute shader):
    osg::Matrix textureToCamView;
    textureToCamView.invert(camViewToTexture);
    cs._texToViewMatrix->set(textureToCamView);

    windDrawable->updateBuffers(cv->getCurrentCamera());

    return cs._stateSet.get();
}
