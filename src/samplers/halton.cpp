
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


// samplers/halton.cpp*
#include "samplers/halton.h"
#include "paramset.h"
#include "rng.h"

namespace pbrt {

// HaltonSampler Local Constants
static PBRT_CONSTEXPR int kMaxResolution = 128;

// HaltonSampler Utility Functions
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y);

// 这个函数的命名完全有问题..MultiplicativeInverse是倒数..这里和倒数没半毛钱关系
// 实际上是 Modular multiplicative inverse

// > 参考：[逆元与取模意义下的除法 - yhang323 - 博客园](https://www.cnblogs.com/wyh323/articles/inverse.html)
//  对于\frac{a}{b} \mod n，如果b在模n下有逆元，那么\frac{a}{b} \mod n = a \cdot b^{-1} \mod n
// 其中  b^{-1} \mod n 是一个标量，可以算出来
// b^{-1} \mod n  =  multiplicativeInverse(b, n);

// 举个例子，比如  a = 18, b = 3, n = 7, \frac{a}{b} mod n = 6
// 这里代入进去， b^{-1} \mod n  =  multiplicativeInverse(3, 7) = 5;
// \frac{18}{3} mod {7} = (18 * 5) mod 7 = 6
// 逆元可以消去除法，转到乘法
static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
    int64_t x, y;

    // ax + ny = gcd(a,n)
    extendedGCD(a, n, &x, &y);
    assert(a * x + n * y == 1); //根据wiki，要使得模拟元存在，gcd必须==1, 也就是a和n必须互素
    return Mod(x, n);
}

// https://zh.wikipedia.org/wiki/%E6%89%A9%E5%B1%95%E6%AC%A7%E5%87%A0%E9%87%8C%E5%BE%97%E7%AE%97%E6%B3%95
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
    if (b == 0) {
        *x = 1;
        *y = 0;
        return;
    }
    int64_t d = a / b, xp, yp;
    extendedGCD(b, a % b, &xp, &yp);
    *x = yp;
    *y = xp - (d * yp);
}

// HaltonSampler Method Definitions
HaltonSampler::HaltonSampler(int samplesPerPixel, const Bounds2i &sampleBounds,
                             bool sampleAtPixelCenter)
    : GlobalSampler(samplesPerPixel), sampleAtPixelCenter(sampleAtPixelCenter) {
    // Generate random digit permutations for Halton sampler

    // 计算扰动表的部分，打表(..
    // 这里就不用看了，不是重点..
    // 具体可以看 https://www.blurredcode.com/2024/01/0f112f8b/#%e6%89%93%e4%b9%b1%e5%ba%8f%e5%88%970%e4%b8%8d%e7%ad%89%e4%ba%8e0%e7%9a%84%e6%83%85%e5%86%b5
    // faure permutation可以打表
    if (radicalInversePermutations.empty()) {
        RNG rng;
        radicalInversePermutations = ComputeRadicalInversePermutations(rng);
    }

    // Find radical inverse base scales and exponents that cover sampling area
    Vector2i res = sampleBounds.pMax - sampleBounds.pMin;
    for (int i = 0; i < 2; ++i) {
        int base = (i == 0) ? 2 : 3;
        int scale = 1, exp = 0;
        while (scale < std::min(res[i], kMaxResolution)) {
            scale *= base;
            ++exp;
        }

        // 假设res = (x,y), 分别计算2^m >= x, 2^n >= y, m,n的值
        // 比如res = {200,150}, kMaxResolution = 128
        // 相当于 {128,128}
        // 那么m=7, n=5, 2^7=128, 3^5=243
        // baseScales = {128,243}
        // baseExps = {7,5}
        baseScales[i] = scale;
        baseExponents[i] = exp;
    }

    // Compute stride in samples for visiting each pixel area
    sampleStride = baseScales[0] * baseScales[1];

    // Compute multiplicative inverses for _baseScales_
    multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
    multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
}

std::vector<uint16_t> HaltonSampler::radicalInversePermutations;
int64_t HaltonSampler::GetIndexForSample(int64_t sampleNum) const {
    if (currentPixel != pixelForOffset) {
        // Compute Halton sample offset for _currentPixel_
        offsetForCurrentPixel = 0;
        if (sampleStride > 1) {
            Point2i pm(Mod(currentPixel[0], kMaxResolution),
                       Mod(currentPixel[1], kMaxResolution));
            for (int i = 0; i < 2; ++i) {
                uint64_t dimOffset =
                    (i == 0)
                        ? InverseRadicalInverse<2>(pm[i], baseExponents[i])
                        : InverseRadicalInverse<3>(pm[i], baseExponents[i]);
                offsetForCurrentPixel +=
                    dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            offsetForCurrentPixel %= sampleStride;
        }
        pixelForOffset = currentPixel;
    }
    return offsetForCurrentPixel + sampleNum * sampleStride;
}

Float HaltonSampler::SampleDimension(int64_t index, int dim) const {
    if (sampleAtPixelCenter && (dim == 0 || dim == 1)) return 0.5f;
    if (dim == 0)
        return RadicalInverse(dim, index >> baseExponents[0]);
    else if (dim == 1)
        return RadicalInverse(dim, index / baseScales[1]);
    else
        return ScrambledRadicalInverse(dim, index,
                                       PermutationForDimension(dim));
}

std::unique_ptr<Sampler> HaltonSampler::Clone(int seed) {
    return std::unique_ptr<Sampler>(new HaltonSampler(*this));
}

HaltonSampler *CreateHaltonSampler(const ParamSet &params,
                                   const Bounds2i &sampleBounds) {
    int nsamp = params.FindOneInt("pixelsamples", 16);
    if (PbrtOptions.quickRender) nsamp = 1;
    bool sampleAtCenter = params.FindOneBool("samplepixelcenter", false);
    return new HaltonSampler(nsamp, sampleBounds, sampleAtCenter);
}

}  // namespace pbrt
