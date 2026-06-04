#ifndef GRADIENT_H_
#define GRADIENT_H_





#include "Vector.hpp"
#include "Mesh.h"
// #include "Field.hpp"


template<typename Tp>
class Field;

template <typename Tp>
class CellField;

template <typename Tp>
class FaceField;

// 定义计算梯度的方式，高斯格林公式、最小二乘
enum class GradientMethod
{
    GAUSS_GREEN,
    LEAST_SQUARES
};




// 计算梯度函数
template <typename Tp>
auto grad(
    const Field<Tp>& field,
    GradientMethod method = GradientMethod::GAUSS_GREEN
) -> CellField<decltype(Tp()* Vector<Scalar>())>
{
    if (!field.isValid())
    {
        std::cerr << "grad() Error: field is not valid." << std::endl;
        throw std::runtime_error("field is not valid.");
    }

    using ULL = unsigned long long;

    using GradType = decltype(Tp()* Vector<Scalar>());  // 梯度类型
    // Field<GradType> resultGradField("Gradient_" + field.getName(), field.getMesh());
    
    CellField<GradType> resultGradField("Gradient_" + field.getName(), field.getMesh());// 取出单元梯度（只给单元赋值）
    resultGradField.initialize();   // 初始化（设置为有效）
    const FaceField<Tp>& faceField = field.getFaceField();
    const std::vector<Cell>& cells = field.getMesh()->getCells();
    const std::vector<Face>& faces = field.getMesh()->getFaces();

    if (method == GradientMethod::GAUSS_GREEN)
    {
        if (field.getMesh()->getDimension() == Mesh::Dimension::TWO_D)
        {
            for (ULL i = 0; i < cells.size(); ++i)
            {
                const Cell& cell = cells[i];
                // 定义总的phi * S_f
                decltype(Tp() * Vector<Scalar>()) total_Phi_Sf{};
                for (ULL j : cell.getFaceIndexes())
                {
                    // 二维需要跳过empty面
                    if (j >= field.getMesh()->getEmptyFaceIndexesPair().first &&
                        j < field.getMesh()->getEmptyFaceIndexesPair().second)
                    {
                        continue;
                    }

                    const Face& face = faces[j];
                    Vector<Scalar> Sf = face.getArea() * face.getNormal();
                    Tp phi = faceField[j];
                    if (face.getOwnerIndex() == i)  // 向外
                    {
                        total_Phi_Sf += phi * Sf;
                    }
                    else
                    {
                        total_Phi_Sf -= phi * Sf;
                    }
                }
                resultGradField[i] = total_Phi_Sf / cell.getVolume();
            }
        }
        else if (field.getMesh()->getDimension() == Mesh::Dimension::THREE_D)
        {
            for (ULL i = 0; i < cells.size(); ++i)
            {
                const Cell& cell = cells[i];
                decltype(Tp() * Vector<Scalar>()) total_Phi_Sf{};
                for (ULL j : cell.getFaceIndexes())
                {
                    const Face& face = faces[j];
                    Vector<Scalar> Sf = face.getArea() * face.getNormal();
                    Tp phi = faceField[j];
                    if (face.getOwnerIndex() == i)  // 向外
                    {
                        total_Phi_Sf += phi * Sf;
                    }
                    else
                    {
                        total_Phi_Sf -= phi * Sf;
                    }
                }
                resultGradField[i] = total_Phi_Sf / cell.getVolume();
            }
        }
        else
        {
            std::cerr << "grad() Error: dimension is not supported." << std::endl;
            throw std::runtime_error("dimension is not supported.");
        }
        return resultGradField;
    }
    else if (method == GradientMethod::LEAST_SQUARES)   // 挖坑
    {
        return resultGradField;
    }
    return resultGradField;
}








/* =============以下放入了Field类的私有函数中============== */


// // 计算场的梯度，默认为高斯格林公式
// template <typename Tp>
// inline auto grad(
//     const Field<Tp>& field,
//     GradientMethod method = GradientMethod::GAUSS_GREEN
// ) -> CellField<decltype(Tp()* Vector<Scalar>())>
// {
//     using ULL = unsigned long long;
//     if (!field.getFaceField().isValid())        // 场必须有效
//     {
//         std::cerr << "Error: field is not valid." << std::endl;
//         throw std::runtime_error("face field is not valid.");
//     }
//     CellField<decltype(Tp()* Vector<Scalar>())> resultField(field.getName(), field.getMesh());
//     resultField.setValue(decltype(Tp() * Vector<Scalar>()){});  // 对齐初始化（设置为有效）
//     // 取出当前面场，遍历每个cell, 计算梯度
//     const FaceField<Tp>& currentFaceField = field.getFaceField();
//     const std::vector<Cell>& cells = field.getMesh()->getCells();
//     const std::vector<Face>& faces = field.getMesh()->getFaces();
//     if (method == GradientMethod::GAUSS_GREEN)
//     {
//         if (field.getMesh()->getDimension() == Mesh::Dimension::TWO_D)
//         {
//             for (int i = 0; i < cells.size(); ++i)  // 计算每个单元的梯度并赋值
//             {
//                 const Cell& cell = cells[i];
//                 // 获取当前单元各个面的id
//                 const std::vector<ULL>& faceIds = cell.getFaceIndexes();
//                 // 定义总的phi * S_f
//                 decltype(Tp() * Vector<Scalar>()) total_Phi_Sf{};
//                 for (ULL j : faceIds)       // 对每个面的Phi * S_f加和
//                 {
//                     // 跳过empty面，如果是else if中三维网格不判断
//                     if (j >= field.getMesh()->getEmptyFaceIndexesPair().first &&
//                         j < field.getMesh()->getEmptyFaceIndexesPair().second)
//                     {
//                         continue;
//                     }
//                     const Face& face = faces[j];
//                     Vector<Scalar> Sf = face.getArea() * face.getNormal();
//                     Tp phi = currentFaceField[j];
//                     total_Phi_Sf += phi * Sf;
//                 }
//                 resultField[i] = total_Phi_Sf / cell.getVolume();
//             }
//         }
//         else if (field.getMesh()->getDimension() == Mesh::Dimension::THREE_D)
//         { 
//             for (int i = 0; i < cells.size(); ++i)  // 计算每个单元的梯度并赋值
//             {
//                 const Cell& cell = cells[i];
//                 // 获取当前单元各个面的id
//                 const std::vector<ULL>& faceIds = cell.getFaceIndexes();
//                 // 定义总的phi * S_f
//                 decltype(Tp() * Vector<Scalar>()) total_Phi_Sf{};
//                 for (ULL j : faceIds)       // 对每个面的Phi * S_f加和
//                 {
//                     const Face& face = faces[j];
//                     Vector<Scalar> Sf = face.getArea() * face.getNormal();
//                     Tp phi = currentFaceField[j];
//                     total_Phi_Sf += phi * Sf;
//                 }
//                 resultField[i] = total_Phi_Sf / cell.getVolume();
//             }
//         }
//         return resultField;
//     }
//     // else if (method == GradientMethod::LEAST_SQUARES)     // 挖坑
//     // {
//     // }
//     return resultField;
// }



#endif // GRADIENT_H_
