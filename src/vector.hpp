#ifndef FILE_VECTOR
#define FILE_VECTOR

#include <iostream>

namespace ASC_bla
{
 
  template <typename T, typename TDIST = std::integral_constant<size_t,1> >
  class VectorView : public VecExpr<VectorView<T,TDIST>>
  {
  protected:
    T * m_data;
    size_t m_size;
    TDIST m_dist;
  public:
    VectorView() = default;
    VectorView(const VectorView &) = default;
    
    template <typename TDIST2>
    VectorView (const VectorView<T,TDIST2> & v2)
      : m_data(v2.data()), m_size(v2.Size()), m_dist(v2.dist()) { }
    
    VectorView (size_t size, T * data)
      : m_data(data), m_size(size) { }
    
    VectorView (size_t size, TDIST dist, T * data)
      : m_data(data), m_size(size), m_dist(dist) { }
    
    template <typename TB>
    VectorView & operator= (const VecExpr<TB> & v2)
    {
      assert (m_size == v2.size());
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = v2(i);
      return *this;
    }

    VectorView & operator= (T scal)
    {
      for (size_t i = 0; i < m_size; i++)
        m_data[m_dist*i] = scal;
      return *this;
    }

    T * data() const { return m_data; }
    size_t size() const { return m_size; }
    auto dist() const { return m_dist; }
    
    T & operator()(size_t i) { return m_data[m_dist*i]; }
    const T & operator()(size_t i) const { return m_data[m_dist*i]; }
    
    auto range(size_t first, size_t next) const {
      assert(first <= next && next <= m_size);
      return VectorView(next-first, m_dist, m_data+first*m_dist);
    }

    auto slice(size_t first, size_t slice) const {
      return VectorView<T,size_t> (m_size/slice, m_dist*slice, m_data+first*m_dist);
    }
      
  };
  

  
  template <typename T>
  class Vector
  {
    size_t size;  // size_t...Datentyp (wie int), aber speziell für Indizes 
    T * data;
    
  public:
    Vector (size_t _size) 
      : size(_size), data(new T[size]) { ; }  //alles nach Doppelpunkt: "size" wird mit _size gesetzt; "data" wird auf T mit Länge size gestellt
    
    Vector (const Vector & v) 
      : Vector(v.Size())      // ✅ ruft den anderen Konstruktor auf!
    {
      *this = v;              // kopiert danach die Daten
    }
    Vector (Vector && v)
      : size(0), data(nullptr)
    {
      std::swap(size, v.size);
      std::swap(data, v.data);
    }

    ~Vector () { delete [] data; }
    
    Vector & operator=(const Vector & v2)
    {
      for (size_t i = 0; i < size; i++)
        data[i] = v2(i);
      return *this;
    }

    Vector & operator= (Vector && v2)
    {
      std::swap(size, v2.size);
      std::swap(data, v2.data);
      return *this;
    }
    
    size_t Size() const { return size; }
    T & operator()(size_t i) { return data[i]; }
    const T & operator()(size_t i) const { return data[i]; }
  };


  /*template <typename ...Args>
  std::ostream & operator<< (std::ostream & ost, const VectorView<Args...> & v)
  template <typename T>
  Vector<T> operator+ (const Vector<T> & a, const Vector<T> & b)
  {
    Vector<T> sum(a.Size());
    for (size_t i = 0; i < a.Size(); i++)
      sum(i) = a(i)+b(i);
    return sum;
  }*/

    // Addition
    template <typename T>
    Vector<T> operator+ (const Vector<T> & a, const Vector<T> & b)
    {
      assert(a.size() == b.size());
      Vector<T> sum(a.size());
      for (size_t i = 0; i < a.size(); i++)
        sum(i) = a(i) + b(i);
      return sum;
    }


template <typename ...Args>
std::ostream & operator<< (std::ostream & ost, const VectorView<Args...> & v)
{
  if (v.size() > 0)
    ost << v(0);
  for (size_t i = 1; i < v.size(); i++)
    ost << ", " << v(i);
  return ost;
}


  template <typename T>
  Vector<T> operator- (const Vector<T> & a, const Vector<T> & b)
  {
      assert(a.size() == b.size());
      Vector<T> minus(a.size());
      for (size_t i = 0; i < a.size(); i++)
        minus(i) = a(i) - b(i);
      return minus;
  }
  
  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const Vector<T> & v)
  {
      if (v.size() > 0)
        ost << v(0);
      for (size_t i = 1; i < v.size(); i++)
        ost << ", " << v(i);
      return ost;
  }
  
}

#endif
