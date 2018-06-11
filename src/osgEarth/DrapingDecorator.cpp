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
#include <osgEarth/DrapingDecorator>
#include <osgEarth/Registry>
#include <osgEarth/Capabilities>
#include <osgEarth/Horizon>
#include <osgEarth/VirtualProgram>
#include <osgEarth/Shaders>
#include <osgEarth/TerrainResources>
#include <osgEarth/ShaderUtils>
#include <osgEarth/LineDrawable>
#include <osgEarth/GLUtils>
#include <osgEarth/Registry>

#include <osg/Texture2D>
#include <osg/Texture2DArray>
#include <osg/BlendFunc>
#include <osg/ShapeDrawable>
#include <osg/AutoTransform>
#include <osg/Depth>
#include <osgUtil/CullVisitor>
#include <osgUtil/LineSegmentIntersector>
#include <osgShadow/ConvexPolyhedron>

#include <stdlib.h> // getenv

#define LC "[DrapingDecorator] "

using namespace osgEarth;

// TODO ITEMS:
// - Address the dangling CameraLocal when a view disappears
// - Clamp maxExt to the frustum far plane extents

DrapingDecorator::DrapingDecorator(const SpatialReference* srs, TerrainResources* resources) :
_unit(-1),
_multisamples(2u),
_maxCascades(4u),
_texSize(1024u),
_mipmapping(false),
_maxHorizonDistance(DBL_MAX),
_debug(false),
_srs(srs),
_resources(resources)
{
    if (::getenv("OSGEARTH_DRAPING_DEBUG"))
        _debug = true;

    const char* c = ::getenv("OSGEARTH_DRAPING_TEXTURE_SIZE");
    if (c)
        setTextureSize((unsigned)atoi(c));

    c = ::getenv("OSGEARTH_DRAPING_MAX_CASCADES");
    if (c)
        setMaxNumCascades((unsigned)atoi(c));

    c = ::getenv("OSGEARTH_DRAPING_MIPMAPPING");
    if (c)
        setUseMipMaps(atoi(c)? true : false);

    c = ::getenv("OSGEARTH_DRAPING_MULTISAMPLES");
    if (c)
        setNumMultiSamples((unsigned)atoi(c));

    c = ::getenv("OSGEARTH_DRAPING_MAX_HORIZON_DISTANCE");
    if (c)
        _maxHorizonDistance = (double)atoi(c);
}

void
DrapingDecorator::setMaxNumCascades(unsigned value)
{
    _maxCascades = osg::clampBetween(value, 1u, 4u);
}

void
DrapingDecorator::setNumMultiSamples(unsigned value)
{
    _multisamples = osg::clampBetween(value, 0u, 4u);
}

void
DrapingDecorator::setTextureSize(unsigned value)
{
    _texSize = osg::clampBetween(value, 256u, 4096u);
}

void
DrapingDecorator::setUseMipMaps(bool value)
{
    _mipmapping = value;
}

void
DrapingDecorator::traverse(osg::NodeVisitor& nv)
{
    bool traversedChildren = false;

    if (nv.getVisitorType() == nv.CULL_VISITOR)
    {
        osgUtil::CullVisitor* cv = dynamic_cast<osgUtil::CullVisitor*>(&nv);
        if (cv)
        {
            const osg::Camera* camera = cv->getCurrentCamera();

            // only proceed if there is geometry to drape.
            // TODO: is this correct? if there's nothing, should be clear out any
            // pre-existing projected texture or set a uniform or something?
            if (_manager.get(camera).getBound().valid())
            {
                // if we don't have a texture unit reserved, do so now.
                if (_unit < 0)
                {
                    reserveTextureImageUnit();
                }

                if (_unit >= 0)
                {
                    // access the draping configuration for this camera:
                    CameraLocal& local = _data.get(camera);

                    // traverse the RTT camera(s) and generate the projective texture
                    local.traverse(cv, *this);

                    // then push the projected texture state and traverse the terrain.
                    cv->pushStateSet(local._terrainSS.get());
                    osg::Group::traverse(nv);
                    cv->popStateSet();

                    traversedChildren = true;
                
                    // debugging
                    if (camera->getName() == "dump")
                        local.dump(camera, *this);

                    local._projMatrixLastFrame = *cv->getProjectionMatrix();
                }
            }
        }
    }

    if (!traversedChildren)
    {
        osg::Group::traverse(nv);
    }
}

void
DrapingDecorator::reserveTextureImageUnit()
{
    if (_unit < 0)
    {
        static Threading::Mutex mutex;
        mutex.lock();

        osg::ref_ptr<TerrainResources> tr;
        if (_unit < 0 && _resources.lock(tr))
        {
            tr->reserveTextureImageUnit(_unit, "Draping");
        }

        mutex.unlock();
    }
}

osg::Node*
DrapingDecorator::getDump()
{
    osg::Node* n = _dump.release();
    _dump = 0L;
    return n;
}

//........................................................................

namespace
{
    /**
     * A camera that will traverse the per-thread DrapingCullSet instead of its own children.
     */
    class DrapingCamera : public osg::Camera
    {
    public:
        DrapingCamera(DrapingManager& dm, const osg::Camera* parentCamera) 
            : osg::Camera(), _parentCamera(parentCamera), _dm(dm)
        {
            setCullingActive( false );
        }

    public: // osg::Node

        void traverse(osg::NodeVisitor& nv)
        {
            DrapingCullSet& cullSet = _dm.get(_parentCamera);
            cullSet.accept( nv );
        }

    protected:
        virtual ~DrapingCamera() { }
        const osg::Camera* _parentCamera;
        DrapingManager& _dm;
    };


    //! Intersect a clip-space ray with a plane.
    bool
    intersectRayWithPlane(const osg::Vec3d& L0,
                          const osg::Vec3d& L1,
                          const osg::Plane& plane,
                          osg::Vec3d& output)
    {
        osg::Vec3d L = L1 - L0;
        L.normalize();

        const osg::Plane::Vec3_type N = plane.getNormal();

        double dist = plane.distance(L0);

#if 0
        osg::Vec3d P0 = L0-N*dist; // point on the plane
        double numer = (P0 - L0) * N; // always -dist? optimize out.
        if (numer == 0)
            return false;
#endif

        double numer = -dist;

        double denom = L * N;
        if (osg::equivalent(denom, 0.0))
            return false;

        double d = numer/denom;

        // enforce "ray" .. fail if intersection is behind L0
        if (d <= 0.0)
            return false;

        output = L0 + L*d;
        return true;
    }

    inline double mix(double a, double b, double t) 
    { 
        return a*(1.0 - t) + b*t;
    }
}

DrapingDecorator::CameraLocal::~CameraLocal()
{
    osg::ref_ptr<osg::Camera> camera;
    if (_token.lock(camera))
    {
        camera->removeObserver(this);
    }
}

void
DrapingDecorator::CameraLocal::objectDeleted(void* ptr)
{
    OE_DEBUG << "CameraLocal::objectDeleted" << std::endl;
    for (unsigned i = 0; i < 4; ++i)
    {
        _cascades[i]._rtt = 0L;
        _terrainSS = 0L;
    }
}

void
DrapingDecorator::CameraLocal::initialize(osg::Camera* camera, DrapingDecorator& decorator)
{
    // set up the auto-delete of orphaned cameras
    _token = camera;
    camera->getOrCreateObserverSet()->addObserver(this);

    unsigned textureWidth;
    unsigned textureHeight;

    unsigned multiSamples;
    osg::Texture::FilterMode minifyFilter;
    osg::Vec4 clearColor;
    bool mipmapping = decorator._mipmapping;

    // if the master cam is a picker, just limit to one cascade with no sampling.
    bool isPickCamera = camera->getName() == "osgEarth::RTTPicker";
    if (isPickCamera)
    {
        textureWidth = osg::minimum(512u, decorator._texSize);
        textureHeight = osg::minimum(512u, decorator._texSize);
        _maxCascades = 2u; // limit to one cascade when picking
        multiSamples = 0u; // no antialiasing allowed
        minifyFilter = osg::Texture::NEAREST; // no texture filtering allowed
        clearColor.set(0,0,0,0);
        mipmapping = false;
    }
    else
    {
        textureWidth = decorator._texSize;
        textureHeight = decorator._texSize;
        _maxCascades = decorator._maxCascades;
        multiSamples = decorator._multisamples;
        minifyFilter = mipmapping? osg::Texture::LINEAR_MIPMAP_LINEAR : osg::Texture::LINEAR;
        clearColor.set(1,1,1,0);
    }

    // Create the shared draping texture.
    osg::Texture2DArray* tex = new osg::Texture2DArray();
    tex->setTextureSize(textureWidth, textureHeight, 4u); //_maxCascades);
    tex->setInternalFormat(GL_RGBA);
    tex->setSourceFormat(GL_RGBA);
    tex->setSourceType(GL_UNSIGNED_BYTE);
    tex->setResizeNonPowerOfTwoHint(false);
    tex->setFilter(tex->MIN_FILTER, minifyFilter);
    tex->setFilter(tex->MAG_FILTER, tex->LINEAR);
    tex->setWrap(tex->WRAP_S, tex->CLAMP_TO_EDGE);
    tex->setWrap(tex->WRAP_T, tex->CLAMP_TO_EDGE);
    tex->setMaxAnisotropy(4.0f);
    
    // set up the global RTT camera state:
    _rttSS = new osg::StateSet();
    osg::StateAttribute::OverrideValue forceOff = osg::StateAttribute::OFF | osg::StateAttribute::PROTECTED | osg::StateAttribute::OVERRIDE;
    osg::StateAttribute::OverrideValue forceOn  = osg::StateAttribute::OFF | osg::StateAttribute::PROTECTED | osg::StateAttribute::OVERRIDE;
    
    // blending:
    _rttSS->setAttributeAndModes(
        new osg::BlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA),
        forceOn);

    // Cannot do this because it will break picking -gw
    //VirtualProgram* rttVP = VirtualProgram::getOrCreate(_rttSS.get());
    //rttVP->setInheritShaders(false);

    for (unsigned i = 0; i < _maxCascades; ++i)
    {
        osg::ref_ptr<osg::Camera> rtt = new DrapingCamera(decorator._manager, camera);

        rtt->setClearColor(clearColor);

        if (decorator._debug && camera->getName() == "dump")
        {
            if (i==0)
                rtt->setClearColor(osg::Vec4(1,0,0,0.15));
            else if (i==1)
                rtt->setClearColor(osg::Vec4(0,1,0,0.15));
            else if (i==2)
                rtt->setClearColor(osg::Vec4(0,0,1,0.15));
            else
                rtt->setClearColor(osg::Vec4(1,0,1,0.15));
        }

        rtt->setGraphicsContext(camera->getGraphicsContext());
        rtt->setReferenceFrame(rtt->ABSOLUTE_RF_INHERIT_VIEWPOINT);
        rtt->setViewport(0, 0, textureWidth, textureHeight);
        rtt->setRenderOrder(rtt->PRE_RENDER);
        rtt->setRenderTargetImplementation(rtt->FRAME_BUFFER_OBJECT);
        rtt->setComputeNearFarMode(rtt->DO_NOT_COMPUTE_NEAR_FAR);
        rtt->setImplicitBufferAttachmentMask(0, 0); // no implicit attachments!
        rtt->setDrawBuffer(GL_FRONT);
        rtt->setReadBuffer(GL_FRONT);

        rtt->attach(
            rtt->COLOR_BUFFER,        // target
            tex,                      // texture to populate
            0u,                       // mipmap level
            i,                        // texture array index
            mipmapping,               // mipmapping
            multiSamples);            // antialiasing multi-samples

        // Note: No depth buffer.

        // Set this so we can detect the RTT camera's parent for the DrapingCamera and
        // for things like auto-scaling, picking, etc.
        // GW: unnecessary unless we want to eventually share RTT cameras?
        rtt->setView(camera->getView());

        // no addChild() because the DrapingCamera will automatically traverse the current cull set

        _cascades[i]._rtt = rtt.get();
    }

    // Set up a stateSet for the terrain that will apply the projected texture.
    _terrainSS = new osg::StateSet();

    // bind the projected texture
    _terrainSS->setTextureAttributeAndModes(decorator._unit, tex, 1);
    _terrainSS->getOrCreateUniform("oe_Draping_tex", osg::Uniform::SAMPLER_2D)->set((int)decorator._unit);
    
    // install the shader program to project a texture on the terrain
    VirtualProgram* drapingShader = VirtualProgram::getOrCreate(_terrainSS.get());
    drapingShader->setName("Draping");
    Shaders shaders;
    shaders.load(drapingShader, shaders.Draping);
}

void
DrapingDecorator::CameraLocal::clear()
{
    //todo - clear out unused RTT cameras/textures
}

void DrapingDecorator::Cascade::expandToInclude(const osg::Vec3d& pointView)
{
    _minX = std::min(_minX, pointView.x());
    _minY = std::min(_minY, pointView.y());
    _maxX = std::max(_maxX, pointView.x());
    _maxY = std::max(_maxY, pointView.y());
}

void DrapingDecorator::Cascade::computeProjection(const osg::Matrix& rttView,
                                                  const osg::Matrix& iCamMVP,
                                                  const osg::Plane& plane,
                                                  double dp, 
                                                  const osg::Vec2d& maxExt,
                                                  double prevMaxY)
{
    // intersect the view frustum's edge vectors with the horizon plane
    osg::Vec3d LL, LR, UL, UR;
    bool LL_ok, LR_ok, UL_ok, UR_ok;
    UL_ok = intersectRayWithPlane(osg::Vec3d(-1, _maxClipY, -1)*iCamMVP, osg::Vec3d(-1, _maxClipY, +1)*iCamMVP, plane, UL);
    UR_ok = intersectRayWithPlane(osg::Vec3d(+1, _maxClipY, -1)*iCamMVP, osg::Vec3d(+1, _maxClipY, +1)*iCamMVP, plane, UR);
    LL_ok = intersectRayWithPlane(osg::Vec3d(-1, _minClipY, -1)*iCamMVP, osg::Vec3d(-1, _minClipY, +1)*iCamMVP, plane, LL);
    LR_ok = intersectRayWithPlane(osg::Vec3d(+1, _minClipY, -1)*iCamMVP, osg::Vec3d(+1, _minClipY, +1)*iCamMVP, plane, LR);

    // next, transform each frustum point into RTT view space (looking down),
    // and expand the orthographic extents to include all valid points            
    _minX = DBL_MAX, _minY = DBL_MAX, _maxX = -DBL_MAX, _maxY = -DBL_MAX;

    // if all frustum edges intersect the horizon plane, calculate the region
    // bounded by the intersected area. Do this by projecting each point into
    // the RTT's view plane.
    if (UL_ok && UR_ok && LL_ok && LR_ok)
    {
        expandToInclude(UL*rttView);
        expandToInclude(UR*rttView);
        expandToInclude(LL*rttView);
        expandToInclude(LR*rttView);
    }

    // if only the LL and LR edges intersect, the view is pitched up over the horizon. 
    // Clamp the Y extent to 0 to prevent RTT-ing behind the camera.
    else if (LL_ok && LR_ok)
    {
        _minY = 0.0;
    }

    // constrain the orthgraphic extents to the visible horizon in all directions.
    if (_minX == DBL_MAX || _minX < -maxExt.x()) _minX = -maxExt.x();
    if (_maxX == -DBL_MAX || _maxX > maxExt.x()) _maxX = maxExt.x();
    if (_minY == DBL_MAX || _minY < -maxExt.y()) _minY = prevMaxY;
    if (_maxY == -DBL_MAX || _maxY > maxExt.y()) _maxY = maxExt.y();

    // If looking forward, clamp the minimum Y to zero so we don't
    // clip geometry in front of us
    if (_minClipY == -1.0 && _minY > 0.0)
        _minY = 0.0;

    else if (_minClipY > -1.0)
        _minY = prevMaxY;

    // finally, compute the projection matrix.
    // TODO: -dp may not always be enough of a near plane...? depends how high up the draped geom is.
    _rttProj.makeOrtho(_minX, _maxX, _minY, _maxY, -dp*4, dp);
}

// Computes the "coverage" of the RTT region in normalized [0..1] clip space.
// If the width and height are 1.0, that means the RTT region will fit exactly 
// within the camera's viewport. For example, a heightNDC of 3.0 means that the
// RTT region is 3x "longer" than the viewport can accommodate; therefore we will need
// multiple cascades.
void DrapingDecorator::Cascade::computeClipCoverage(const osg::Matrix& rttView, const osg::Matrix& camMVP)
{
    // inverse of the RTT MVP matrix:
    osg::Matrix rttMVP = rttView * _rttProj;
    osg::Matrix rttMVPInv;
    rttMVPInv.invert(rttMVP);

    // matrix to transform from RTT clip to camera clip:
    osg::Matrix rttClipToCamClip = rttMVPInv * camMVP;

    osg::Vec3d winLL = osg::Vec3d(-1, _minClipY, +1) * rttClipToCamClip;
    osg::Vec3d winLR = osg::Vec3d(+1, _minClipY, +1) * rttClipToCamClip;
    osg::Vec3d winUL = osg::Vec3d(-1, _maxClipY, +1) * rttClipToCamClip;
    osg::Vec3d winUR = osg::Vec3d(+1, _maxClipY, +1) * rttClipToCamClip;

    // width and height [0..1]
    _widthNDC = std::max(winUR.x() - winUL.x(), winLR.x() - winLL.x())*0.5 + 0.5;
    _heightNDC = (std::max(winUL.y(), winUR.y()) - std::min(winLL.y(), winLR.y()))*0.5 + 0.5;
}

#define MAXABS4(A,B,C,D) \
    osg::maximum(fabs(A), osg::maximum(fabs(B), osg::maximum(fabs(C),fabs(D))))

void DrapingDecorator::CameraLocal::constrainMaxExtToFrustum(const osg::Matrix& iCamMVP, 
                                                             const osg::Matrix& rttView, 
                                                             osg::Vec2d& maxExt)
{
    // transform camera clip space into RTT view space (the maxExt space)
    osg::Matrix camProjToRttView = iCamMVP * rttView;

    const bool constrainY = false;

    if (constrainY)
    {
        // contrain both dimensions - in practice, this might not be necessary
        // because the far clip plane is always (?) beyond the horizon plane due to 
        // bounding-sphere culling. Also, doing this may cause the cascade sizes to
        // "jump" which might be an undesirable visual effect. -gw
        osg::Vec3d farLL = osg::Vec3d(-1,-1,+1)*camProjToRttView;
        osg::Vec3d farLR = osg::Vec3d(+1,-1,+1)*camProjToRttView;
        osg::Vec3d farUL = osg::Vec3d(-1,+1,+1)*camProjToRttView;
        osg::Vec3d farUR = osg::Vec3d(+1,+1,+1)*camProjToRttView;
        maxExt.x() = std::min(maxExt.x(), MAXABS4(farLL.x(), farLR.x(), farUL.x(), farUR.x()));
        maxExt.y() = std::min(maxExt.y(), MAXABS4(farLL.y(), farLR.y(), farUL.y(), farUR.y()));
    }
    else
    {
        osg::Vec3d farLL = osg::Vec3d(-1,-1,+1)*camProjToRttView;
        osg::Vec3d farUR = osg::Vec3d(+1,+1,+1)*camProjToRttView;
        maxExt.x() = osg::minimum(maxExt.x(), osg::maximum(fabs(farLL.x()), fabs(farUR.x())));
    }
}

void rttY_to_camClipY(double rttY, double& outClipY)
{
    osg::Matrix rttViewToCamClip;
    // = invRttView * camMVP;
}

void
DrapingDecorator::CameraLocal::traverse(osgUtil::CullVisitor* cv, DrapingDecorator& decorator)
{    
    osg::Camera* camera = cv->getCurrentCamera();

    // first time through, intiailize the RTT cameras.
    if (_rttSS.valid() == false)
    {
        initialize(camera, decorator);
    }

    // establish the camera view vectors:
    const osg::Matrix& camView = camera->getViewMatrix();
    osg::Vec3d camEye, camCenter, camUp;
    camView.getLookAt(camEye, camCenter, camUp);
    osg::Vec3d camLook = camCenter-camEye;
    camLook.normalize();

    // camera modelview matrix (world -> view)
    const osg::Matrix& camMV = *cv->getModelViewMatrix();

    // camera projection matrix (from previous frame)
    const osg::Matrix& camProj = _projMatrixLastFrame;

    // camera world -> clip
    osg::Matrixd camMVP = camMV * camProj;

    // camera clip -> world
    osg::Matrixd iCamMVP;
    iCamMVP.invert(camMVP);

    // horizon plane (world space) - the plane passing through the
    // ellipsoid's visible horizon in all directions from the camera.
    osg::Plane horizonPlane;

    // distance to the visible horizon (in any direction)
    double dh;

    // shortest distance from the camera to the horizon plane (straight down)
    double dp;
    
    // Maximum theorectical extent of the draping region; i.e. the distance
    // from the camera to the visible horizon projected onto the horizon plane.
    // This the largest possible extent we will need for RTT.
    osg::Vec2d maxExt;

    if (decorator._srs->isGeographic())
    {
        Horizon* horizon = Horizon::get(*cv);
        horizon->getPlane(horizonPlane);
        dh = horizon->getDistanceToVisibleHorizon();
        dp = horizonPlane.distance(camEye);

#if 0
        if (dh > decorator._maxHorizonDistance)
        {
            dh = decorator._maxHorizonDistance;
        }

        // intersect the terrain at the bottom of the view frustum.
        double minY = -maxExt.y();
        osg::Vec3d lookedAt = osg::Vec3d(0, -1, +1) * iCamMVP; // far plane, bottom center of view.
        osg::ref_ptr<osgUtil::LineSegmentIntersector> lsi = new osgUtil::LineSegmentIntersector(camEye, lookedAt);
        osgUtil::IntersectionVisitor iv(lsi.get());
        decorator.accept(iv);
        if (lsi->containsIntersections())
        {
            osg::Vec3d terrainPt = lsi->getFirstIntersection().getWorldIntersectPoint();
            osg::ref_ptr<Horizon> h = new Horizon();
            h->setEllipsoid(osg::EllipsoidModel(terrainPt.length(), terrainPt.length()));
            h->getPlane(horizonPlane);
            dh = h->getDistanceToVisibleHorizon();            
            dp = horizonPlane.distance(camEye);
        }
#endif
    }
    else
    {
        // in projected mode, the horizon is at an infinite distance, so
        // we need to simulate it.
        dp = osg::maximum(camEye.z(), 100.0);
        dh = sqrt(2.0*6356752.3142*dp + dp*dp);
        horizonPlane.set(osg::Vec3d(0,0,1), dp);
    }

    // project visible horizon distance into the horizon plane:
    double m = sqrt(dh*dh - dp*dp);    
    maxExt.set(m, m);

    // Create a view matrix that looks straight down at the horizon plane form the eyepoint.
    // This will be our view matrix for all RTT draping cameras.
    osg::Matrix rttView;
    osg::Vec3d rttLook = -horizonPlane.getNormal();
    osg::Vec3d camLeft = camUp ^ camLook;
    osg::Vec3d rttUp = rttLook ^ camLeft;
    rttView.makeLookAt(camEye, camEye + rttLook, rttUp);

#if 1
    // intersect the terrain at the bottom of the view frustum.
    double y0 = -maxExt.y();
    double y1 =  maxExt.y();
    osg::Vec3d lookedAt = osg::Vec3d(0, -1, +1) * iCamMVP; // far plane, bottom center of view.
    osg::ref_ptr<osgUtil::LineSegmentIntersector> lsi = new osgUtil::LineSegmentIntersector(camEye, lookedAt);
    osgUtil::IntersectionVisitor iv(lsi.get());
    decorator.accept(iv);
    if (lsi->containsIntersections())
    {
        osg::Vec3d terrainPt = lsi->getFirstIntersection().getWorldIntersectPoint();
        terrainPt = terrainPt * rttView;
        y0 = terrainPt.y();
    }
#endif

    // so far, maxExt is a theoretical max. Now we can constrain it based on
    // the actual far clip plane.
    constrainMaxExtToFrustum(iCamMVP, rttView, maxExt);

    // camera view -> world
    osg::Matrix iCamMV;
    iCamMV.invert(camMV);

    // xform from clip [-1..1] to texture [0..1] space
    static const osg::Matrix clipToTex =
        osg::Matrix::translate(1.0, 1.0, 1.0) *
        osg::Matrix::scale(0.5, 0.5, 0.5);

    // Start by computing the full extent of the RTT region. We will use the results
    // of this math to decide how many cascades we need to use. (We are using 
    // cascade[0] here as a temporary workspace, but will re-use it later for the actual
    // highest resolution cascade).
    _cascades[0]._minClipY = -1.0;
    _cascades[0]._maxClipY =  1.0;
    _cascades[0].computeProjection(rttView, iCamMVP, horizonPlane, dp, maxExt, -maxExt.y());

    // Next compute the extent, in pixels, of the full RTT. If it's larger than our
    // texture cascade size, we may need multiple cascases.
    _cascades[0].computeClipCoverage(rttView, camMVP);

    // Prepare to write to the texture matrix array.
    // (Use the gloal maxCascades here regardless of this CameraLocal's maxCascades)
    ArrayUniform texMat("oe_Draping_texMatrix", osg::Uniform::FLOAT_MAT4, _terrainSS.get(), decorator._maxCascades);

    {
        // RTT view to Camera Clip
        osg::Matrix iRttView;
        iRttView.invert(rttView);
        osg::Matrix rttViewToCamClip = iRttView * camMVP;

        _numCascades = 4;
        float f;
        osg::Uniform* u = dynamic_cast<osg::Uniform*>(Registry::instance()->dataStore().fetch(Registry::instance(), "ff"));
        if (u) u->get(f); else f = 0.25;

        double sy = y1-y0;
        int c = 3;
        _cascades[c]._maxY = y1;
        _cascades[c]._minY = y0+sy*f;
        _cascades[c]._minClipY = (osg::Vec3d(0, _cascades[c]._minY, -dp) * rttViewToCamClip).y();
        _cascades[c]._maxClipY = (osg::Vec3d(0, _cascades[c]._maxY, -dp) * rttViewToCamClip).y();
        
        --c;
        sy = _cascades[c+1]._minY - y0;
        _cascades[c]._maxY = _cascades[c+1]._minY;
        _cascades[c]._minY = y0+sy*f;
        _cascades[c]._minClipY = (osg::Vec3d(0, _cascades[c]._minY, -dp) * rttViewToCamClip).y();
        _cascades[c]._maxClipY = (osg::Vec3d(0, _cascades[c]._maxY, -dp) * rttViewToCamClip).y();
        
        --c;
        sy = _cascades[c+1]._minY - y0;
        _cascades[c]._maxY = _cascades[c+1]._minY;
        _cascades[c]._minY = y0+sy*f;
        _cascades[c]._minClipY = (osg::Vec3d(0, _cascades[c]._minY, -dp) * rttViewToCamClip).y();
        _cascades[c]._maxClipY = (osg::Vec3d(0, _cascades[c]._maxY, -dp) * rttViewToCamClip).y();

        --c;
        _cascades[c]._maxY = _cascades[c+1]._minY;
        _cascades[c]._minY = y0;
        _cascades[c]._minClipY = (osg::Vec3d(0, _cascades[c]._minY, -dp) * rttViewToCamClip).y();
        _cascades[c]._maxClipY = (osg::Vec3d(0, _cascades[c]._maxY, -dp) * rttViewToCamClip).y();
    }

#if 0
    // HARD-CODED CASCADE limits. We determined these hueristically. Later if necessary
    // we can attempt to compute optimal values using the results of the computeClipCoverage
    // method.
    {
        // heightNDC is the Y-extent of the draping region's normalized [0..1] clip space.
        // For example, if heightHDC is 3.0, that means the draping region is 3x deeper 
        // than the texture size we're using, so we need multiple cascades to cover it.
        // The following values we determined hueristacally.
        double h = _cascades[0]._heightNDC;

        _cascades[0]._minClipY = -1.00;
        _cascades[0]._maxClipY = -0.95; //-1.00 + 1.0/h;

        _cascades[1]._minClipY = _cascades[0]._maxClipY;
        _cascades[1]._maxClipY = -0.15;

        _cascades[2]._minClipY = _cascades[1]._maxClipY;
        _cascades[2]._maxClipY =  0.35;

        _cascades[3]._minClipY = _cascades[2]._maxClipY;
        _cascades[3]._maxClipY =  1.00;

        _numCascades =
            h >= 4.5 && _maxCascades >= 4 ? 4 :
            h >= 3.0 && _maxCascades >= 3 ? 3 :
            h >= 1.5 && _maxCascades >= 2 ? 2 : 1;

        // For the "camera looking mostly down" case, split the frustum in half to get
        // maximum coverage. Another hueristic tweak.
        double camDotRtt = camLook * rttLook;
        if (_numCascades == 2 && camDotRtt > 0.5)
        {
            double m = _cascades[0]._maxClipY;
            _cascades[0]._maxClipY = _cascades[1]._minClipY = mix(m, 0, camDotRtt);
        }

        _cascades[_numCascades-1]._maxClipY = 1.0;
    }
#endif

    // For each active cascade, configure its RTT camera and build the
    // texture projection matrix
    unsigned i;
    for (i = 0; i < _numCascades; ++i)
    {
        Cascade& cascade = _cascades[i];
        osg::Camera* rtt = cascade._rtt.get();

        // Only do this we have more than one cascade, because we already computed
        // the "overall" region earlier :)
#if 0
        if (_numCascades > 1u)
        {
            double prevMaxY = i == 0 ? -maxExt.y() : _cascades[i-1]._maxY;
            cascade.computeProjection(rttView, iCamMVP, horizonPlane, dp, maxExt, prevMaxY);
        }
#else
        double prevMaxY = i == 0 ? y0 : _cascades[i-1]._maxY;
        cascade.computeProjection(rttView, iCamMVP, horizonPlane, dp, maxExt, prevMaxY);
#endif

        // configure the RTT camera's matrices:
        rtt->setViewMatrix(rttView);
        rtt->setProjectionMatrix(cascade._rttProj);

        // Create the texture matrix that will transform the RTT frame into texture [0..1] space.
        // Doing this on the CPU avoids precision errors on the GPU.
        texMat.setElement(i, iCamMV * rttView * cascade._rttProj * clipToTex);
    }

    if (i < 4)
    {
        // install a "marker" matrix that tells the shader we're past the final cascade.
        static osg::Matrix marker(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
        texMat.setElement(i, marker);
    }

    // traverse and write to the texture.
    if (_numCascades > 0u)
    {
        cv->pushStateSet(_rttSS.get());
        for (unsigned i = 0; i < _numCascades; ++i)
        {
            _cascades[i]._rtt->accept(*cv);
        }
        cv->popStateSet(); // _rttSS
    }
}

void
DrapingDecorator::CameraLocal::dump(const osg::Camera* cam, DrapingDecorator& decorator)
{
    static const char* fn = "DrapingDecoratorDump.osgb";
 
    decorator._dump = new osg::Group();

    // Main camera:
    {
        osgShadow::ConvexPolyhedron ph;
        ph.setToUnitFrustum();
        osg::Matrix mvp = cam->getViewMatrix() * cam->getProjectionMatrix();
        osg::Matrix imvp; imvp.invert(mvp);
        ph.transform( imvp, mvp );
        ph.dumpGeometry(0,0,0,fn);
        osg::ref_ptr<osg::Node> camNode = osgDB::readRefNodeFile(fn);
        camNode->setName("camera");
        decorator._dump->addChild(camNode.get());
    }

    // RTT cameras
    for (unsigned i=0; i<_numCascades; ++i)
    {
        const osg::Camera* rtt = _cascades[i]._rtt.get();
        osgShadow::ConvexPolyhedron ph;
        ph.setToUnitFrustum();
        osg::Matrix mvp = rtt->getViewMatrix() * rtt->getProjectionMatrix();
        osg::Matrix imvp; imvp.invert(mvp);
        ph.transform( imvp, mvp );
        ph.dumpGeometry(0, 0, 0, fn, osg::Vec4(1, 1, 0, 1), osg::Vec4(1, 1, 0, 0.25));
        osg::ref_ptr<osg::Node> camNode = osgDB::readRefNodeFile(fn);
        camNode->setName("rtt");
        decorator._dump->addChild(camNode.get());
    }
}
