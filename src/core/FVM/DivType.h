#ifndef DIVTYPE_H_
#define DIVTYPE_H_



namespace fvm {
    /**
     * @brief 对流项迎风格式
     */
    enum class DivType {
        FUD,    // 一阶迎风
        SUD,    // 二阶迎风
        CD,     // 中心差分
        MINMOD,
        MUSCL
    };
}


#endif // DIVTYPE_H_
