#include <sstream>
#include <pybind11/pybind11.h>

#include "vector.hpp"

using namespace ASC_bla;
namespace py = pybind11;

PYBIND11_MODULE(bla, m) {
    m.doc() = "Basic linear algebra module"; // optional module docstring
    
    py::class_<Vector<double>> (m, "Vector")
      .def(py::init<size_t>(),
           py::arg("size"), "create vector of given size")
      .def("__len__", &Vector<double>::size,
           "return size of vector")
      
      .def("__setitem__", [](Vector<double> & self, int i, double v) {
        if (i < 0) i += self.size();
        if (i < 0 || i >= self.size()) throw py::index_error("vector index out of range");
        self(i) = v;
      })
      .def("__getitem__", [](Vector<double> & self, int i) { return self(i); })
      
      .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
      {
        size_t start, stop, step, n;
        if (!inds.compute(self.size(), &start, &stop, &step, &n))
          throw py::error_already_set();
        self.range(start, stop).slice(0,step) = val;
      })
      
      .def("__add__", [](Vector<double> & self, Vector<double> & other)
      { return Vector<double> (self+other); })

      .def("__rmul__", [](Vector<double> & self, double scal)
      { return Vector<double> (scal*self); })
      
      .def("__str__", [](const Vector<double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      })

     .def(py::pickle(
        [](Vector<double> & self) { // __getstate__
            //return a tuple that fully encodes the state of the object 
          return py::make_tuple(self.size(),
                                py::bytes((char*)(void*)&self(0), self.size()*sizeof(double)));
        },
        [](py::tuple t) { // __setstate__
          if (t.size() != 2)
            throw std::runtime_error("should be a 2-tuple!");

          Vector<double> v(t[0].cast<size_t>());
          py::bytes mem = t[1].cast<py::bytes>();
          std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.size()*sizeof(double));
          return v;
        }))
    ;
}


/*
#include <sstream>
#include <cstring>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>          
#include <pybind11/functional.h>

#include "vector.hpp"
#include "matrix.hpp"    //enthält auch matrixexpr.hpp

namespace py = pybind11;
using namespace ASC_bla;

static inline size_t norm_index(long i, size_t n) {
    if (i < 0) i += static_cast<long>(n);
    if (i < 0 || static_cast<size_t>(i) >= n)
        throw py::index_error("index out of range");
    return static_cast<size_t>(i);
}

PYBIND11_MODULE(bla, m) {
    m.doc() = "Basic linear algebra module (vectors & matrices)";

    //Vector<double> 
    py::class_<Vector<double>>(m, "Vector")
        .def(py::init<size_t>(), py::arg("size"), "create vector of given size")
        .def("__len__", &Vector<double>::size, "return size of vector")

        .def("__setitem__", [](Vector<double>& self, long i, double v) {
            auto ii = norm_index(i, self.size());
            self(ii) = v;
        })
        .def("__getitem__", [](Vector<double>& self, long i) {
            auto ii = norm_index(i, self.size());
            return self(ii);
        })

        .def("__setitem__", [](Vector<double>& self, py::slice inds, double val) {
            size_t start, stop, step, n;
            if (!inds.compute(self.size(), &start, &stop, &step, &n))
                throw py::error_already_set();
            self.range(start, stop).slice(0, step) = val;
        })

        .def("__add__", [](const Vector<double>& a, const Vector<double>& b) {
            if (a.size() != b.size()) throw std::runtime_error("vector sizes differ");
            return Vector<double>(a + b);
        })
        .def("__rmul__", [](const Vector<double>& v, double s) {
            return Vector<double>(s * v);
        })

        .def("__str__", [](const Vector<double>& self) {
            std::stringstream ss; ss << self; return ss.str();
        })

        .def(py::pickle(
            [](const Vector<double>& self) { //__getstate__
                return py::make_tuple(self.size(),
                    py::bytes(reinterpret_cast<const char*>(&self(0)),
                              self.size() * sizeof(double)));
            },
            [](py::tuple t) { //__setstate__
                if (t.size() != 2) throw std::runtime_error("bad state for Vector");
                Vector<double> v(t[0].cast<size_t>());
                py::bytes mem = t[1].cast<py::bytes>();
                std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                            v.size() * sizeof(double));
                return v;
            }));

    //Matrix<double> Binding
    py::class_<Matrix<double>>(m, "Matrix",
        R"doc(Matrix of doubles with simple expression templates (add/mul) and Gauss-Jordan inverse.)doc")
        //Konstruktoren
        .def(py::init<size_t, size_t>(), py::arg("rows"), py::arg("cols"),
             "Create empty matrix with given shape")

        //Hilfskonstruktor aus verschachtelten Listen
        .def_static("from_list", [](const std::vector<std::vector<double>>& rows) {
            if (rows.empty()) throw std::runtime_error("from_list: empty data");
            const size_t r = rows.size();
            const size_t c = rows[0].size();
            for (const auto& row : rows)
                if (row.size() != c) throw std::runtime_error("from_list: ragged rows");
            Matrix<double> M(r, c);
            for (size_t i = 0; i < r; ++i)
                for (size_t j = 0; j < c; ++j)
                    M(i, j) = rows[i][j];
            return M;
        }, "Create Matrix from a list of lists")

        //Basic Info
        .def("rows", &Matrix<double>::Rows)
        .def("cols", &Matrix<double>::Cols)
        .def_property_readonly("shape", [](const Matrix<double>& M) {
            return py::make_tuple(M.Rows(), M.Cols());
        })

        //2D Indexing: M[i, j] / M[i, j] = x
        .def("__getitem__", [](Matrix<double>& self, py::tuple ij) {
            if (ij.size() != 2) throw py::index_error("matrix indexing needs 2 indices");
            auto i = norm_index(ij[0].cast<long>(), self.Rows());
            auto j = norm_index(ij[1].cast<long>(), self.Cols());
            return self(i, j);
        })
        .def("__setitem__", [](Matrix<double>& self, py::tuple ij, double v) {
            if (ij.size() != 2) throw py::index_error("matrix indexing needs 2 indices");
            auto i = norm_index(ij[0].cast<long>(), self.Rows());
            auto j = norm_index(ij[1].cast<long>(), self.Cols());
            self(i, j) = v;
        })

        //Addition: A + B 
        .def("__add__", [](const Matrix<double>& A, const Matrix<double>& B) {
            if (A.Rows() != B.Rows() || A.Cols() != B.Cols())
                throw std::runtime_error("matrix shapes differ in +");
            return Matrix<double>(A + B);   
        }, py::is_operator())

        //Multiplikation: A * B
        .def("__mul__", [](const Matrix<double>& A, const Matrix<double>& B) {
            if (A.Cols() != B.Rows())
                throw std::runtime_error("shapes not aligned for *");
            return Matrix<double>(A * B);   
        }, py::is_operator())

        //Inverse (nur quadratisch)
        .def("inverse", [](const Matrix<double>& A) {
            if (A.Rows() != A.Cols())
                throw std::runtime_error("inverse: matrix must be square");
            return Inverse(A);
        }, "Return A^{-1} using Gauss–Jordan with partial pivoting")

        //nach Python-Listen exportieren 
        .def("tolist", [](const Matrix<double>& M) {
            std::vector<std::vector<double>> out(M.Rows(), std::vector<double>(M.Cols()));
            for (size_t i = 0; i < M.Rows(); ++i)
                for (size_t j = 0; j < M.Cols(); ++j)
                    out[i][j] = M(i, j);
            return out;
        })

        .def("__str__", [](const Matrix<double>& M) {
            std::stringstream ss; ss << M; return ss.str();
        })

        .def(py::pickle(
            [](const Matrix<double>& self) { // __getstate__
                const size_t r = self.Rows(), c = self.Cols();
                std::string blob; blob.resize(r * c * sizeof(double));
                std::memcpy(blob.data(), &self(0, 0), blob.size());
                return py::make_tuple(r, c, py::bytes(blob));
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 3) throw std::runtime_error("bad state for Matrix");
                const size_t r = t[0].cast<size_t>();
                const size_t c = t[1].cast<size_t>();
                Matrix<double> M(r, c);
                py::bytes mem = t[2].cast<py::bytes>();
                std::memcpy(&M(0, 0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                            r * c * sizeof(double));
                return M;
            })
        );

    //Optional - einfach nur bequem
    m.def("matmul", [](const Matrix<double>& A, const Matrix<double>& B) {
        if (A.Cols() != B.Rows()) throw std::runtime_error("shapes not aligned");
        return Matrix<double>(A * B);
    }, "Matrix product A*B");

    m.def("matadd", [](const Matrix<double>& A, const Matrix<double>& B) {
        if (A.Rows() != B.Rows() || A.Cols() != B.Cols())
            throw std::runtime_error("matrix shapes differ");
        return Matrix<double>(A + B);
    }, "Matrix sum A+B");
}


*/