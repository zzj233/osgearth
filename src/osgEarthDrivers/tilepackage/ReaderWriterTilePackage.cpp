/* -*-c++-*- */
/* osgEarth - Dynamic map generation toolkit for OpenSceneGraph
* Copyright 2016 Pelican Mapping
* http://osgearth.org
*
* osgEarth is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
* IN THE SOFTWARE.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*/
#include "TilePackageOptions"
#include "BundleReader"

#include <osgEarth/TileSource>
#include <osgEarth/Registry>
#include <osgEarth/URI>

#include <osg/Notify>
#include <osgDB/FileNameUtils>
#include <osgDB/FileUtils>
#include <osgDB/Registry>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>

#define BUNDLE_PACKET_SIZE 128

using namespace osgEarth;
using namespace osgEarth::Drivers;

#define LC "[ReaderWriterTilePackage] "

std::string padLeft(std::string value, unsigned int length)
{    
    std::stringstream ss;
    if (value.size() < length)
    {        
        for (unsigned int i = 0; i < (length - value.size()); i++)
        {
            ss << "0";
        }
        ss << value;
        return ss.str();
    }
    else
    {
        return value;
    }        
}

class TilePackageSource : public TileSource
{
public:
    TilePackageSource( const TileSourceOptions& options ) :
      TileSource( options ),
      _options( options ),
      _profileConf( ProfileOptions() )
    {        
    }

    // override
    Status initialize( const osgDB::Options* dbOptions )
    {
        URI url = _options.url().value();

		OE_NOTICE << LC << "Initial URL: " << url.full() << std::endl;

        _dbOptions = Registry::instance()->cloneOrCreateOptions( dbOptions );        

        // establish a profile if we don't already have one:
        if ( !getProfile() )
        {
            const Profile* profile = NULL;

            if ( _profileConf.isSet() )
            {
                profile = Profile::create( _profileConf.get() );
            }            
            else
            {
                // finally, fall back on lat/long
                profile = osgEarth::Registry::instance()->getGlobalGeodeticProfile();
            }
            setProfile( profile );
        }

        BundleReader reader("D:/geodata/simdis/arcgis_tile_package/Midnight_Canvas_4326/CanvasMaps_Midnight/Layers/_alllayers/L08/R0080C0080.bundle");        

        return STATUS_OK;
    }

    // override
    int getPixelsPerTile() const
    {
        // TODO:  Get size from conf
        return 256;
    }

    // override
    osg::Image* createImage(const TileKey& key, ProgressCallback* progress)
    {
        // Try to figure out which bundle file the incoming tilekey is in.
        unsigned int numWide, numHigh;
        getProfile()->getNumTiles(key.getLevelOfDetail(), numWide, numHigh);

        std::stringstream buf;
        buf << _options.url()->full() << "/_alllayers/";
        buf << "L" << padLeft(toString<unsigned int>(key.getLevelOfDetail()), 2) << "/";
        
        unsigned int colOffset = floor(numWide / 128);
        unsigned int rowOffset = floor(numHigh / 128);

        OE_NOTICE << "Offsets=" << colOffset << "x" << rowOffset << std::endl;        

        buf << "R" << padLeft(toHex(rowOffset), 4) << "C" << padLeft(toHex(colOffset), 4);
        buf << ".bundle";


        std::string bundleFile = buf.str();
        OE_NOTICE << "Bundle file=" << bundleFile << std::endl;

        if (osgDB::fileExists(bundleFile))
        {
            BundleReader reader(bundleFile);
            osg::Image* result = reader.readImage(key);
            /*
            if (result)
            {
                std::stringstream outFile;
                outFile << key.getLevelOfDetail() << "_" << key.getTileX() << "_" << key.getTileY() << ".png";
                osgDB::writeImageFile(*result, outFile.str());
            }
            */
            return result;
        }
        else
        {
            OE_NOTICE << bundleFile << " doesn't exist" << std::endl;
        }
        
        return 0;


        //BundleReader reader("D:/geodata/simdis/arcgis_tile_package/Midnight_Canvas_4326/CanvasMaps_Midnight/Layers/_alllayers/L00/R0000C0000.bundle");
        //return reader.readImage(0);
    }

    // override
    virtual std::string getExtension() const 
    {
        // TODO:  Get from config
        return "png";
    }

private:
    const TilePackageOptions _options;
    optional<ProfileOptions> _profileConf;    
    osg::ref_ptr<osgDB::Options> _dbOptions;
};


class TilePackageTileSourceFactory : public TileSourceDriver
{
public:
    TilePackageTileSourceFactory()
    {
        supportsExtension( "osgearth_tilepackage", "TilePackage" );
    }

    virtual const char* className() const
    {
        return "TilePackage ReaderWriter";
    }

    virtual ReadResult readObject(const std::string& file_name, const Options* options) const
    {
        if ( !acceptsExtension(osgDB::getLowerCaseFileExtension( file_name )))
            return ReadResult::FILE_NOT_HANDLED;

        return new TilePackageSource( getTileSourceOptions(options) );
    }
};

REGISTER_OSGPLUGIN(osgearth_tilepackage, TilePackageTileSourceFactory)


