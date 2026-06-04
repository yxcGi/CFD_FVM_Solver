#ifndef DIVERGENCE_H_
#define DIVERGENCE_H_

#include "Vector.hpp"
#include "Mesh.h"

template<typename Tp>
class Field;

template <typename Tp>
class CellField;

template <typename Tp>
class FaceField;

// 散度计算
template <typename Tp>
auto divergence(const FaceField<Tp>& field) -> CellField<decltype(Tp()& Vector<Scalar>())> {
    if (!field.isValid()) {
        std::cerr << "divergence(const Field<Tp>& field): field is invalid." << std::endl;
        throw std::invalid_argument("field is invalid.");
    }
    using ULL = unsigned long long;
    using DivType = decltype(Tp()& Vector<Scalar>());

    CellField<DivType> resultDivField("Divergence_" + field.getName(), field.getMesh());
    resultDivField.setValue(DivType{});
    const std::vector<Cell>& cells = field.getMesh()->getCells();
    const std::vector<Face>& faces = field.getMesh()->getFaces();


    if (field.getMesh()->getDimension() == Mesh::Dimension::TWO_D) {
        for (ULL i = 0; i < cells.size(); ++i) {
            const Cell& cell = cells[i];
            DivType div{};
            for (ULL j : cell.getFaceIndexes()) {
                // 二维需要跳过empty面
                if (j >= field.getMesh()->getEmptyFaceIndexesPair().first &&
                    j < field.getMesh()->getEmptyFaceIndexesPair().second)
                {
                    continue;
                }

                const Face& face = faces[j];
                Vector<Scalar> Sf = face.getArea() * face.getNormal();
                Tp phi = field[j];
                if (face.getOwnerIndex() == i) {
                    div += phi & Sf;
                }
                else {
                    div -= phi & Sf;
                }
            }
            resultDivField[i] = div;
        }
    }
    else if (field.getMesh()->getDimension() == Mesh::Dimension::THREE_D) {
        for (ULL i = 0; i < cells.size(); ++i) {
            const Cell& cell = cells[i];
            DivType div{};
            for (ULL j : cell.getFaceIndexes()) {
                const Face& face = faces[j];
                Vector<Scalar> Sf = face.getArea() * face.getNormal();
                Tp phi = field[j];
                if (face.getOwnerIndex() == i) {
                    div += phi & Sf;
                }
                else {
                    div -= phi & Sf;
                }
            }
            resultDivField[i] = div;
        }
    }
    return resultDivField;
}

#endif // DIVERGENCE_H_