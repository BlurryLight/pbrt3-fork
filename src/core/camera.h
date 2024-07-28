
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "film.h"

namespace pbrt {

// Camera Declarations
class Camera {
  public:
    // Camera Interface
    Camera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
           Float shutterClose, Film *film, const Medium *medium);
    virtual ~Camera();
    virtual Float GenerateRay(const CameraSample &sample, Ray *ray) const = 0;
    virtual Float GenerateRayDifferential(const CameraSample &sample,
                                          RayDifferential *rd) const;
    virtual Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
    virtual void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
    virtual Spectrum Sample_Wi(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf, Point2f *pRaster,
                               VisibilityTester *vis) const;

    // Camera Public Data
    AnimatedTransform CameraToWorld;
    const Float shutterOpen, shutterClose;
    Film *film;
    const Medium *medium;
};

struct CameraSample {
    Point2f pFilm;
    Point2f pLens;
    Float time;
};

inline std::ostream &operator<<(std::ostream &os, const CameraSample &cs) {
    os << "[ pFilm: " << cs.pFilm << " , pLens: " << cs.pLens <<
        StringPrintf(", time %f ]", cs.time);
    return os;
}

class ProjectiveCamera : public Camera {
  public:
    // ProjectiveCamera Public Methods
    ProjectiveCamera(const AnimatedTransform &CameraToWorld,
                     const Transform &CameraToScreen,
                     const Bounds2f &screenWindow, Float shutterOpen,
                     Float shutterClose, Float lensr, Float focald, Film *film,
                     const Medium *medium)
        : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
          CameraToScreen(CameraToScreen) {
        // Initialize depth of field parameters
        lensRadius = lensr;
        focalDistance = focald;

        // Compute projective camera transformations

        // Compute projective camera screen transformations
        // read from bottom to upper
        // first translate
        // then scale to (1/..  ,经过前两步缩放到了(0,1)范围，也就到了NDC
        // then scale the file->fullResolution

        // ScreenSpace: 经过了投影变化，坐标系还是和ViewSpace是一样的
        // 注意，ScreenSpace和NDCSpace的Y轴是相反的，所以第二步的scale是(1,-1,1)
        // ScreenSpace -> NDCSpace
        // (pMin.xy) --> (0,1)  (pMin.x, pMax.y) -> (0,0)
        // (pMax.xy) --> (1,0)  (pMax.x, pMin.y) -> (1,1)

        //  ![Untitled-1-2024-07-28-19-25-44](https://img.blurredcode.com/img/Untitled-1-2024-07-28-19-25-44.png?x-oss-process=style/compress)
        // 首先把左上角对齐到(0,0)点，需要平移 (-pMin.x, -pMax.y)
        // 然后Y轴需要反转，X,Y需要缩放到(0,1),
        // 所以X轴的缩放为 (pMax.x - pMin.x)
        // Y轴的缩放为(pMax.y - pMin.y) * (-1) //反转Y轴
        ScreenToRaster =
            Scale(film->fullResolution.x, film->fullResolution.y, 1) *
            Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x),
                  1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
            Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));

        // nori和pbrt的一个关键差别点是
        // 对于ScreenToCamera矩阵，
        // nori以width为1，height为 1 / aspect
        // pbrt以height为1，width为 1 * aspect
    
       // pbrt的1 / (screenWindow.pMax.x - screenWindow.pMin.x) 约等于 (0.5 / aspect,0.5,1)
       // 而nori的这部分是Eigen::DiagonalMatrix<float, 3>(Vector3f(-0.5f, -0.5f * aspect, 1.0f))
       // Translation部分也要做部分变换
        float aspect = (float)film->fullResolution.x / film->fullResolution.y;
        auto noriScreenWindow =
            Bounds2f(screenWindow.pMin / aspect, screenWindow.pMax / aspect);
        auto noriScreenToRaster =
            Scale(film->fullResolution.x, film->fullResolution.y, 1) *
            Scale(1 / (noriScreenWindow.pMax.x - noriScreenWindow.pMin.x),
                  1 / (noriScreenWindow.pMin.y - noriScreenWindow.pMax.y), 1)
                   *
            Translate(Vector3f(-noriScreenWindow.pMin.x, -noriScreenWindow.pMax.y, 0));

        // RasterToScreen = Inverse(noriScreenToRaster);

        RasterToScreen = Inverse(ScreenToRaster);
        RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;
    }

  protected:
    // ProjectiveCamera Protected Data
    Transform CameraToScreen, RasterToCamera;
    Transform ScreenToRaster, RasterToScreen;
    Float lensRadius, focalDistance;
};

}  // namespace pbrt

#endif  // PBRT_CORE_CAMERA_H
