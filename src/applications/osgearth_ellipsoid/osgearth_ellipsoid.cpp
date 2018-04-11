#include <osgViewer/Viewer>
#include <osg/CoordinateSystemNode>
#include <osgViewer/ViewerEventHandlers>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osgDB/ReadFile>
#include <osg/CullFace>
#include <osgEarth/VirtualProgram>


const char* VS =
"#version " GLSL_VERSION_STR "\n"

"uniform mat4  osg_ViewMatrixInverse;\n"
"uniform mat4  osg_ViewMatrix;\n"

"out vec3 world;\n"

"void main(void)\n"
"{\n"
"   gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * gl_Vertex;\n"
"   mat4 modelMatrix = osg_ViewMatrixInverse * gl_ModelViewMatrix;\n"
"   vec4 vert = modelMatrix  * gl_Vertex;\n"
"   world = vert.xyz;\n"
"    gl_FrontColor = gl_Color;\n"
"}\n";

const char* FS =
"#version " GLSL_VERSION_STR "\n"

"uniform vec3 oe_eye; \n"
"uniform float radius; \n"
"in vec3 world; \n"

"float raySphereIntersect(vec3 r0, vec3 rd, vec3 s0, float sr) {\n"
"float a = dot(rd, rd); \n"
"   vec3 s0_r0 = r0 - s0; \n"
"   float b = 2.0 * dot(rd, s0_r0); \n"
"   float c = dot(s0_r0, s0_r0) - (sr * sr); \n"
"   if (b*b - 4.0*a*c < 0.0) {\n"
"        return -1.0; \n"
"    }\n"
"    return (-b - sqrt((b*b) - 4.0*a*c)) / (2.0*a); \n"
"}\n"
"void main(void)\n"
"{\n"
"    vec3 direction = normalize(world - oe_eye);\n"
"    float hit = raySphereIntersect(oe_eye, direction, vec3(0, 0, 0), radius);\n"
"    if (hit <= 0.0) discard;\n"
"    gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);\n"
"}\n";


osg::Node* makeEllipsoid(double radius)
{
    osg::Geometry* geometry = new osg::Geometry;
    
    osg::Vec3Array* verts = new osg::Vec3Array;
    geometry->setVertexArray(verts);

    double minX = -radius;
    double minY = -radius;
    double minZ = -radius;
    double maxX = radius;
    double maxY = radius;
    double maxZ = radius;

    osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1, 0, 0, 1));
    
    geometry->setColorArray(colors, osg::Array::BIND_OVERALL);

    verts->push_back(osg::Vec3(minX, minY, minZ));
    verts->push_back(osg::Vec3(maxX, minY, minZ));
    verts->push_back(osg::Vec3(minX, minY, maxZ));
    verts->push_back(osg::Vec3(maxX, minY, maxZ));
    
    verts->push_back(osg::Vec3(minX, maxY, minZ));
    verts->push_back(osg::Vec3(maxX, maxY, minZ));
    verts->push_back(osg::Vec3(minX, maxY, maxZ));
    verts->push_back(osg::Vec3(maxX, maxY, maxZ));

    osg::DrawElementsUShort* de = new osg::DrawElementsUShort(GL_TRIANGLES);
    //front
    de->push_back(0); de->push_back(1); de->push_back(2);
    de->push_back(1); de->push_back(3); de->push_back(2);
    
    // back
    de->push_back(5); de->push_back(4); de->push_back(7);
    de->push_back(4); de->push_back(6); de->push_back(7);

    // right
    de->push_back(1); de->push_back(5); de->push_back(3);
    de->push_back(5); de->push_back(7); de->push_back(3);

    // left
    de->push_back(4); de->push_back(0); de->push_back(6);
    de->push_back(6); de->push_back(0); de->push_back(2);

    // bottom
    de->push_back(4); de->push_back(1); de->push_back(0);
    de->push_back(4); de->push_back(5); de->push_back(1);

    // top
    de->push_back(2); de->push_back(3); de->push_back(6);
    de->push_back(6); de->push_back(3); de->push_back(7);

    geometry->addPrimitiveSet(de);

    osg::Program* program = new osg::Program;
    //program->addShader(osgDB::readShaderFile(osg::Shader::VERTEX, "ellipsoid_vs.glsl"));
    //program->addShader(osgDB::readShaderFile(osg::Shader::FRAGMENT, "ellipsoid_fs.glsl"));
    program->addShader(new osg::Shader(osg::Shader::VERTEX, VS));
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, FS));

    geometry->getOrCreateStateSet()->setAttributeAndModes(new osg::CullFace(osg::CullFace::FRONT));

    geometry->getOrCreateStateSet()->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF);

    geometry->getOrCreateStateSet()->getOrCreateUniform("radius", osg::Uniform::FLOAT)->set((float)radius);

    geometry->getOrCreateStateSet()->setAttributeAndModes(program);

    return geometry;
}

class CameraCullCallback : public osg::NodeCallback
{
    virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
    {
        osg::Camera* camera = static_cast<osg::Camera*>(node);
        osg::Vec3d eye, center, up;
        camera->getViewMatrixAsLookAt(eye, center, up);
        camera->getOrCreateStateSet()->getOrCreateUniform("oe_eye", osg::Uniform::FLOAT_VEC3)->set(osg::Vec3(eye));
        traverse(node, nv);
    }
};

int
main(int argc, char** argv)
{
    osg::ArgumentParser arguments(&argc,argv);

    // set up a viewer:
    osgViewer::Viewer viewer(arguments);

    osg::Group* root = new osg::Group;
    root->addChild(makeEllipsoid(osg::WGS_84_RADIUS_EQUATOR));
    viewer.setSceneData(root);

    // add some stock OSG handlers:
    viewer.addEventHandler(new osgViewer::StatsHandler());
    viewer.addEventHandler(new osgViewer::WindowSizeHandler());
    viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));

    viewer.setCameraManipulator(new osgGA::TrackballManipulator());

    viewer.getCamera()->setCullCallback(new CameraCullCallback());

    viewer.realize();

    while (!viewer.done())
    {       
        viewer.frame();
    }

    return 0;
}
