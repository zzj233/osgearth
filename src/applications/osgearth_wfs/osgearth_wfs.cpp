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

#include <osgEarth/XmlUtils>
#include <osgEarth/TileKey>
#include <osgEarth/Registry>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgDB/FileNameUtils>
#include <osg/Math>

using namespace osgEarth;

#define INDEX_HEADER_SIZE 16
#define INDEX_SIZE 5

unsigned int  computeOffset(const std::vector<char>& buffer) {
    unsigned int sum = 0;
    for (unsigned int i = 0; i < buffer.size(); i++) {
        char v = buffer[i];
        //std::cout << (unsigned int)v << std::endl;
        sum += ((unsigned int)v & 0xff) * pow(2, 8 * i);
    }
    return sum;
}

class BundleReader
{
public:
    BundleReader(const std::string& directory):
        _directory(directory),
        _bundleDim(128) // This is in the conf but I think is the same for all of them.
    {    
    }

    /**
     * Reads the index of a buffer file.
     */
    void readIndex(const std::string& filename, std::vector<int>& index)
    {        
        std::ifstream input(filename.c_str(), std::ifstream::binary);
        char header[INDEX_HEADER_SIZE];
        input.read(header, INDEX_HEADER_SIZE);
        while (input.good()) {
            std::vector<char> buffer;
            buffer.resize(5);
            if (input.read(&buffer[0], INDEX_SIZE))
            {
                int offset = computeOffset(buffer);
                index.push_back(offset);
            }
        }
    }

    void readTile(std::istream& buffer, int offset, std::ostream& output)
    {        
        OSG_NOTICE << "Offset " << offset << std::endl;
        buffer.seekg(offset, std::ios::beg);
        std::vector<char> sizeBuffer;
        sizeBuffer.resize(4);
        buffer.read(&sizeBuffer[0], 4);
        int size = computeOffset(sizeBuffer);
        if (size > 0)
        {
            OSG_NOTICE << "Size=" << size << std::endl;

            char* image = new char[size];
            buffer.read(image, size);
            output.write(image, size);
            delete[] image;
        }
        else
        {
            OSG_NOTICE << "Size is 0" << std::endl;
        }
    }

    unsigned int hexFromString(const std::string& input)
    {
        unsigned int result;
        std::stringstream ss;
        ss << std::hex << input;
        ss >> result;
        return result;
    }

    osg::Image* getImage(TileKey &key)
    {
        std::string base = "D:/geodata/simdis/arcgis_tile_package/Midnight_Canvas_4326/CanvasMaps_Midnight/Layers/_alllayers/L08/R0080C0180";
        //std::string base = "D:/geodata/simdis/arcgis_tile_package/Midnight_Canvas_4326/CanvasMaps_Midnight/Layers/_alllayers/L08/RffffCffff";
        std::string bundleFile = base + ".bundle";
        std::string indexFile = base + ".bundlx";
        
        readIndex(indexFile, _index);
        OSG_NOTICE << "Read " << _index.size() << " indices" << std::endl;

        std::string baseName = osgDB::getNameLessExtension(osgDB::getSimpleFileName(bundleFile));
        std::cout << "Basename=" << baseName << std::endl;

        /*
        unsigned int colOffset, rowOffset;       
        std::stringstream ss;
        ss << std::hex << baseName.substr(1, 4);
        std::cout << "Row=" << ss.str();
        ss >> rowOffset;

        ss.str("");
        ss.clear();
        ss << std::hex << baseName.substr(6, 4);
        std::cout << "Col=" << ss.str();
        ss >> colOffset;
        */
        unsigned int rowOffset = hexFromString(baseName.substr(1, 4));
        unsigned int colOffset = hexFromString(baseName.substr(6, 4));

        std::cout << "row=" << rowOffset << ", col=" << colOffset << std::endl;



        std::ifstream bundle(bundleFile.c_str(), std::ofstream::binary);

        for (unsigned int i = 0; i < _index.size(); i++)
        {
            unsigned int row = floor(float(i) / _bundleDim);                
            unsigned int x = colOffset + row;
            unsigned int y = rowOffset + i - (row * _bundleDim);

            std::cout << "Writing index " << i << " (" << x << ", " << y << ")" << std::endl;
            std::stringstream outName;
            outName << i << ".png";

            std::stringstream buf;
            readTile(bundle, _index[i], buf);
            if (buf.str().size() > 0)
            {
                std::string data = buf.str();
                std::ofstream out(outName.str().c_str(), std::ofstream::binary);
                out.write(data.c_str(), data.size());
            }
        }

        return 0;       
    }

    int _bundleDim;

    std::string _directory;
    std::vector< int > _index;
    int _indexSize;
};

int
main(int argc, char** argv)
{
    BundleReader reader("D:/geodata/simdis/arcgis_tile_package/Midnight_Canvas_4326/CanvasMaps_Midnight/Layers");

    const Profile* profile = osgEarth::Registry::instance()->getGlobalGeodeticProfile();
    
    osg::Image* image = reader.getImage(TileKey(0, 0, 0, profile));
    if (image) {
        osgDB::writeImageFile(*image, "image.png");
    }
}
