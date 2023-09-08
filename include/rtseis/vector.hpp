#ifndef RTSEIS_VECTOR_HPP
#define RTSEIS_VECTOR_HPP
#include <vector>
#include <memory>
#include <boost/align.hpp>
namespace RTSeis
{
template<class T = double>
/// @brief Defines a vector for use in RTSeis.  This a bit more general
///        than standard C++ vectors and provides simple functions.
class Vector
{
private:
    using DataTypeT = std::vector<T, boost::alignment::aligned_allocator<T, 64>>;
public:
    using iterator = typename DataTypeT::iterator;
    using const_iterator = typename DataTypeT::const_iterator;
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Vector();
    /// @brief Copy constructor.
    Vector(const Vector &v);
    /// @brief Move constructor.
    Vector(Vector &&v) noexcept;
    /// @brief Constructs from a given vector.
    explicit Vector(const std::vector<T> &vector); 
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    Vector& operator=(const Vector &v);
    /// @brief Move assignment.
    Vector& operator=(Vector &&v) noexcept;
    /// @}

    /// @brief True indicates that the vector is empty.
    [[nodiscard]] bool empty() const noexcept;
    /// @brief Resizes the vector to the given size.
    void resize(size_t n);
    /// @brief Resizes the vector to teh given size and fills with the value.
    void resize(size_t n, T value);
    /// @brief Reserves space for the vector.
    void reserve(size_t n);
    /// @result The size of the vector.
    [[nodiscard]] size_t size() const noexcept;
    /// @result A pointer to the internal memory array.
    [[nodiscard]] T *data() noexcept;
    /// @result A pointer to the internal memory array.
    [[nodiscard]] const T *data() const noexcept;

    static int getAlignment() noexcept;

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;

    [[nodiscard]] T& operator[](size_t index);
    [[nodiscard]] T& operator[](size_t index) const;
    [[nodiscard]] T& at(size_t index);
    [[nodiscard]] T& at(size_t index) const;
 
    /// @name Destructors
    /// @{

    /// @brief Releases memory.
    void clear() noexcept;
    /// @brief Destructor. 
    ~Vector();
    /// @}
private:
    class VectorImpl;
    std::unique_ptr<VectorImpl> pImpl;
};
}
#endif
