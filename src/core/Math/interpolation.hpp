#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_





#include <iostream>
#include <type_traits>
#include "Vector.hpp"
#include "Tensor.hpp"


namespace interpolation
{
    enum class Scheme
    {
        // 线性差值
        LINEAR,
        // 算数平均
        ARITHMETIC_MEAN,
        // 调和平均
        HARMONIC_MEAN,
    };
}


namespace util {

    /**
     * @brief 插值类函数对象，用于两点间的差值
     * @tparam T 数据类型
     */
    template <typename Tp>
    class Interpolation;



    // 模板特例化，仅支持三种数据类型的插值：Scalar, Vector<Scalar>, Tensor<Scalar>
    template <>
    class Interpolation<Scalar>
    {
        using Scalar = double;
    public:
        // 传两个值就是算数平均
        // T operator()(const T& phi1, const T& phi2);
        Scalar operator()(
            const Scalar phi1, const Scalar phi2,
            interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
            Scalar alpha = 0.5) const;
    };

    template <>
    class Interpolation<Vector<Scalar>>
    {
        using Scalar = double;
    public:
        // 传两个值就是算数平均
        // T operator()(const T& phi1, const T& phi2);
        Vector<Scalar> operator()(
            const Vector<Scalar>& phi1, const Vector<Scalar>& phi2,
            interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
            Scalar alpha = 0.5) const;
    };

    template <>
    class Interpolation<Tensor<Scalar>>
    {
        using Scalar = double;
    public:
        // 传两个值就是算数平均
        // T operator()(const T& phi1, const T& phi2);
        Tensor<Scalar> operator()(
            const Tensor<Scalar>& phi1, const Tensor<Scalar>& phi2,
            interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
            Scalar alpha = 0.5) const;
    };


    /**
 * @brief 标量插值函数接口。
 *
 * @param phi1 第一个位置的标量值。
 * @param phi2 第二个位置的标量值。
 * @param scheme 插值格式。
 * @param alpha 插值系数。默认 alpha = 0.5。
 *
 * @return 插值后的标量值。
 */
    Scalar interpolate(
        Scalar phi1,
        Scalar phi2,
        interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
        Scalar alpha = 0.5);


    /**
     * @brief 向量插值函数接口。
     *
     * @param phi1 第一个位置的向量值。
     * @param phi2 第二个位置的向量值。
     * @param scheme 插值格式。
     * @param alpha 插值系数。默认 alpha = 0.5。
     *
     * @return 插值后的向量值。
     */
    Vector<Scalar> interpolate(
        const Vector<Scalar>& phi1,
        const Vector<Scalar>& phi2,
        interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
        Scalar alpha = 0.5);


    /**
     * @brief 张量插值函数接口。
     *
     * @param phi1 第一个位置的张量值。
     * @param phi2 第二个位置的张量值。
     * @param scheme 插值格式。
     * @param alpha 插值系数。默认 alpha = 0.5。
     *
     * @return 插值后的张量值。
     */
    Tensor<Scalar> interpolate(
        const Tensor<Scalar>& phi1,
        const Tensor<Scalar>& phi2,
        interpolation::Scheme scheme = interpolation::Scheme::LINEAR,
        Scalar alpha = 0.5);






    inline double Interpolation<double>::operator()(const Scalar phi1, const Scalar phi2, interpolation::Scheme scheme, Scalar alpha) const
    {
        // phi1 * (1 - alpha) + phi2 * alpha
        if (scheme == interpolation::Scheme::LINEAR)
        {
            return phi1 + alpha * (phi2 - phi1);
        }
        else if (scheme == interpolation::Scheme::HARMONIC_MEAN)
        {
            return phi1 * phi2 / (phi2 * alpha + phi1 * (1 - alpha));
        }
        else if (scheme == interpolation::Scheme::ARITHMETIC_MEAN)
        {
            return (phi1 + phi2) * 0.5;
        }
        std::cerr << "interpolation::Scheme not supported!" << std::endl;
        throw std::runtime_error("interpolation::Scheme not supported!");
    }

    inline Vector<Scalar> Interpolation<Vector<Scalar>>::operator()(const Vector<Scalar>& phi1, const Vector<Scalar>& phi2, interpolation::Scheme scheme, Scalar alpha) const
    {
        if (scheme == interpolation::Scheme::LINEAR)
        {
            return (1 - alpha) * phi1 + alpha * phi2;
        }
        else if (scheme == interpolation::Scheme::ARITHMETIC_MEAN)
        {
            return (phi1 + phi2) * 0.5;
        }
        std::cerr << "interpolation::Scheme not supported!" << std::endl;
        throw std::runtime_error("interpolation::Scheme not supported!");
    }


    inline Tensor<Scalar> Interpolation<Tensor<Scalar>>::operator()(const Tensor<Scalar>& phi1, const Tensor<Scalar>& phi2, interpolation::Scheme scheme, Scalar alpha) const
    {
        if (scheme == interpolation::Scheme::LINEAR)
        {
            return (1 - alpha) * phi1 + alpha * phi2;
        }
        else if (scheme == interpolation::Scheme::ARITHMETIC_MEAN)
        {
            return (phi1 + phi2) * 0.5;
        }
        std::cerr << "interpolation::Scheme not supported!" << std::endl;
        throw std::runtime_error("interpolation::Scheme not supported!");
    }


    inline Scalar interpolate(
        const Scalar phi1,
        const Scalar phi2,
        const interpolation::Scheme scheme,
        const Scalar alpha)
    {
        return Interpolation<Scalar>{}(phi1, phi2, scheme, alpha);
    }


    inline Vector<Scalar> interpolate(
        const Vector<Scalar>& phi1,
        const Vector<Scalar>& phi2,
        const interpolation::Scheme scheme,
        const Scalar alpha)
    {
        return Interpolation<Vector<Scalar>>{}(phi1, phi2, scheme, alpha);
    }


    inline Tensor<Scalar> interpolate(
        const Tensor<Scalar>& phi1,
        const Tensor<Scalar>& phi2,
        const interpolation::Scheme scheme,
        const Scalar alpha)
    {
        return Interpolation<Tensor<Scalar>>{}(phi1, phi2, scheme, alpha);
    }


}

#endif // INTERPOLATION_H_
