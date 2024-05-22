#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <future>
template <typename T>
class Matrix {
private:
    int lines;
    int columns;
    T** matrix;
public:
    T convert_to(const std::string &str)
    {
        std::stringstream ss;
        T num;
        ss <<str;
        ss >> std::noskipws >> num;
        return num;
    }
    void Print()
    {
        for (int i = 0; i < lines; ++i)
        {
            std::cout << "[";
            for (int j = 0; j < columns; ++j)
            {
                if (j != columns - 1) { std::cout << matrix[i][j] << " "; }
                else { std::cout << matrix[i][j]; }
            }
            std::cout << "]\n";
        }
        std::cout << "\n";
    };
    void Init()
    {
        std::cout << "Fill out the matrix!\n";
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                std::cin >> matrix[i][j];
            }
        }
    };
    void WriteToAFile()
    {
        std::string FileName;
        std::cout << "Enter the name of file:";
        std::cin >> FileName;
        std::string endFileName="/home/durininkirous/Homework/"+FileName;
        std::ofstream out;
        out.open(endFileName);
        for (int i = 0; i < lines; ++i) {
            out <<"[";
            for (int j = 0; j < columns; ++j) {
                out <<matrix[i][j];
                if (j != columns - 1) { out <<" "; }
            }
            out <<"]\n";
        }
        out.close();
    }
    void FirstElementaryTransformation(int FirstLine, int SecondLine)
    {
        T tmp;
        for (int i = 0; i < columns; ++i)
        {
            tmp = matrix[FirstLine][i];
            matrix[FirstLine][i] = matrix[SecondLine][i];
            matrix[SecondLine][i] = tmp;
        }
    }
    template <typename T1>
    void SecondElementaryTransformation(int ChangeLine, T1 MultiPlier)
    {
        for (int i = 0; i < columns; ++i)
        {
            matrix[ChangeLine][i] *= MultiPlier;
        }
    }
    template <typename T1>
    void ThirdElementaryTransformation(int TakeLine, T1 MultiPlier, int ChangeLine)
    {
        for (int i = 0; i < columns; ++i)
        {
            matrix[ChangeLine][i] += matrix[TakeLine][i] * MultiPlier;
        }
    }
    static Matrix InitializationByZeros(int lines, int columns)
    {
        return Matrix{lines, columns};
    }
    static Matrix InitializationByUnits(int lines, int columns)
    {
        if (lines!=columns) {throw "Error! Identity matrix can only be square";}
        else
        {
            Matrix NewMatrix{lines, columns};
            for (int i=0;i<NewMatrix.lines;++i)
            {
                for (int j=0;j<NewMatrix.columns;++j)
                {
                    if (i==j) {NewMatrix.matrix[i][j]=1;}
                }
            }
            return NewMatrix;
        }
    }
    T Determinant()
    {
        T det = 0;
        if (lines != columns)
        {
            std::cout << "There is no determinant for such matrix\n";
            return 0;
        }
        else {
            if (lines==1)
            {
                return matrix[0][0];
            }
            for (int i = 0; i < columns; i++)
            {
                Matrix<T> NewMatrix(*this, 0, i);
                det += matrix[0][i] * pow(-1, i) * NewMatrix.Determinant();
            }
            return det;
        }
    }
    void Transposition()
    {
        T tmp;
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                tmp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = tmp;
            }
        }
    }
    int GetLines()
    {
        return lines;
    }
    int GetColumns()
    {
        return columns;
    }
    T** GetMatrix()
    {
        return matrix;
    }
    Matrix()
    {
        std::cout << "Enter the number of lines:"; std::cin >> lines;
        std::cout << "Enter  the number of columns:"; std::cin >> columns;
        matrix = new T* [lines];
        for (int i = 0; i < lines; ++i)
        {
            matrix[i] = new T[columns];
        }
        Init();
    };
    Matrix(int lines, int columns, T** OtherMatrix) //Собрать матрицу из уже существующих элементов
    {
        this->lines = lines;
        this->columns = columns;
        T** matrix = new T*[lines];
        for (int i = 0; i < lines; ++i)
        {
            matrix[i] = new T[columns];
        }
        for (int i=0;i<lines;i++)
        {
            for (int j = 0; j < columns; ++j)
            {
                matrix[i][j] = OtherMatrix[i][j];
            }
        }
        this->matrix = matrix;
    }
    Matrix(int lines, int columns) //Собрать нулевую матрицу из lines*columns элементов
    {
        this->lines = lines;
        this->columns = columns;
        T** matrix = new T* [lines];
        for (int i = 0; i < lines; ++i)
        {
            matrix[i] = new T[columns];
        }
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                matrix[i][j] = 0;
            }
        }
        this->matrix = matrix;
    }
    Matrix(const Matrix& p)
    {
        lines = p.lines;
        columns = p.columns;
        T** matrix = new T*[lines];
        for (int i = 0; i < lines; i++)
        {
            matrix[i] = new T[columns];
        }
        for (int i=0;i<lines;++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                matrix[i][j] = p.matrix[i][j];
            }
        }
        this->matrix = matrix;
    }
    Matrix(const Matrix p, int line, int column) //Собрать матрицу на основе переданной, но в которой вычеркнуты строка line и столбец column
    {
        int k = 0;
        int m = 0;
        this->lines = p.lines - 1;
        this->columns = p.columns - 1;
        T** matrix = new T* [lines];
        for (int i = 0; i < lines; ++i)
        {
            matrix[i] = new T[columns];
        }
        for (int i = 0; i < p.lines; ++i)
        {
            m = 0;
            if (i!= line) {
                for (int j = 0; j < p.columns; ++j)
                {
                    if (j!= column) {
                        matrix[k][m] = p.matrix[i][j];
                        m += 1;
                    }
                }
                k += 1;
            }
        }
        this->matrix = matrix;
    }
    Matrix(std::string FileName)//Собрать матрицу из файла, нужно передать имя файла (Уже существует файл MatrixExample.txt, но можно создать свою)
    {
        std::fstream in;
        in.open(FileName);
        std::string line;
        std::string garb="[ ;]";
        int lines = 0, columns = 0;
        while (std::getline(in,line))
        {
            ++lines;
            std::stringstream ss(line);
            std::string token;
            columns=0;
            while (ss >> token)
            {
                    ++columns;
            }
        }
        T** matrix = new T* [lines];
        for (int i=0;i<lines;++i)
        {
            matrix[i] = new T[columns];
        }
        int CurrLine = 0;
        int CurrColumn=0;
        in.close();
        in.open(FileName);
        while (std::getline(in,line))
        {
            CurrColumn = 0;
            line.erase(line.find('['),line.find('[')+1);
            line.erase(line.find(']'),line.find(']'));
            std::stringstream ss(line);
            std::string token;
            ss.str(line);
            while (ss >> token)
            {
                matrix[CurrLine][CurrColumn]=(convert_to(token));
                ++CurrColumn;
            }
            ++CurrLine;
            ss.clear();
        }
        this->lines = lines;
        this->columns = columns;
        this->matrix=matrix;
        in.close();
    };
    ~Matrix()
    {
        for (int i = 0; i < lines; ++i)
        {
            delete[] matrix[i];
        }
        delete[] matrix;
    }
    template <typename T1>
    auto operator + (Matrix<T1> Second)-> Matrix<decltype(matrix[0][0]+Second.GetMatrix()[0][0])>
    {
        using MatrixSum = Matrix<decltype(matrix[0][0]+Second.GetMatrix()[0][0])>;
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            MatrixSum ReturnMatrix(lines, columns);
            for (int i = 0; i < lines; ++i)
            {
                for (int j = 0; j < columns; ++j)
                {
                    ReturnMatrix.GetMatrix()[i][j]=matrix[i][j] + Second.GetMatrix()[i][j];
                }
            }
            return ReturnMatrix;
        }
        else {
            std::cerr << "Error! The matrices must be of the same dimension. Returning a zero matrix with the same dimension as the first one"<<std::endl;
            return MatrixSum{ lines,columns };
        }
    }
    template <typename T1>
    auto Parallel_sum(Matrix<T1> Second) -> Matrix<decltype(matrix[0][0] + Second.GetMatrix()[0][0])>
    {
        using MatrixSum = Matrix<decltype(matrix[0][0] + Second.GetMatrix()[0][0])>;
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            MatrixSum ReturnMatrix(lines, columns);
            std::vector<std::thread> threads;
            for (int i = 0; i < lines; ++i)
            {
                threads.emplace_back([i, &ReturnMatrix, this, &Second]() 
                {
                    for (int j = 0; j < columns; ++j)
                    {
                            ReturnMatrix.GetMatrix()[i][j] = matrix[i][j] + Second.GetMatrix()[i][j];
                    }
                });
            }
            for (auto& thread : threads)
            {
                thread.join();
            }
            return ReturnMatrix;
        }
        else {
            std::cerr << "Error! The matrices must be of the same dimension. Returning a zero matrix with the same dimension as the first one" << std::endl;
            return MatrixSum{ lines, columns };
        }
    }
    template <typename T1>
    auto Parallel_sum(Matrix<T1> Second, int blockSize) -> Matrix<decltype(matrix[0][0] + Second.GetMatrix()[0][0])>
    {
        using MatrixSum = Matrix<decltype(matrix[0][0] + Second.GetMatrix()[0][0])>;
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            MatrixSum ReturnMatrix(lines, columns);
            std::vector<std::future<void>> futures;
            for (int i = 0; i < lines; i += blockSize)
            {
                futures.emplace_back(std::async(std::launch::async, [i, blockSize, &ReturnMatrix, this, &Second]()
                    {
                        int end = std::min(i + blockSize, lines);
                        for (int row = i; row < end; ++row)
                        {
                            for (int j = 0; j < columns; ++j)
                            {
                                ReturnMatrix.GetMatrix()[row][j] = matrix[row][j] + Second.GetMatrix()[row][j];
                            }
                        }
                    }));
            }

            for (auto& future : futures)
            {
                future.get();
            }
            return ReturnMatrix;
        }
        else {
            std::cerr << "Error! The matrices must be of the same dimension. Returning a zero matrix with the same dimension as the first one" << std::endl;
            return MatrixSum{ lines, columns };
        }
    }
    template <typename T1>
    auto operator - (Matrix<T1> Second)-> Matrix<decltype(matrix[0][0]-Second.GetMatrix()[0][0])>
    {
        using MatrixSub = Matrix<decltype(matrix[0][0]-Second.GetMatrix()[0][0])>;
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            MatrixSub ReturnMatrix(lines, columns);
            for (int i = 0; i < lines; ++i)
            {
                for (int j = 0; j < columns; ++j)
                {
                    ReturnMatrix.GetMatrix()[i][j]=matrix[i][j] - Second.GetMatrix()[i][j];
                }
            }
            return ReturnMatrix;
        }
        else {
            std::cerr << "Error! The matrices must be of the same dimension. Returning a zero matrix with the same dimension as the first one"<<std::endl;
            return MatrixSub{ lines,columns };
        }
    }
    template <typename T1>
    auto Parallel_razn(Matrix<T1> Second) -> Matrix<decltype(matrix[0][0] - Second.GetMatrix()[0][0])>
    {
        using MatrixSum = Matrix<decltype(matrix[0][0] - Second.GetMatrix()[0][0])>;
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            MatrixSum ReturnMatrix(lines, columns);
            std::vector<std::thread> threads;
            for (int i = 0; i < lines; ++i)
            {
                threads.emplace_back([i, &ReturnMatrix, this, &Second]()
                    {
                        for (int j = 0; j < columns; ++j)
                        {
                            ReturnMatrix.GetMatrix()[i][j] = matrix[i][j] - Second.GetMatrix()[i][j];
                        }
                    });
            }
            for (auto& thread : threads)
            {
                thread.join();
            }
            return ReturnMatrix;
        }
        else {
            std::cerr << "Error! The matrices must be of the same dimension. Returning a zero matrix with the same dimension as the first one" << std::endl;
            return MatrixSum{ lines, columns };
        }
    }
    template <typename T1>
    auto Parallel_razn(Matrix<T1> Second, int blockSize) -> Matrix<decltype(matrix[0][0] - Second.GetMatrix()[0][0])>
    {
        using MatrixSum = Matrix<decltype(matrix[0][0] - Second.GetMatrix()[0][0])>;
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            MatrixSum ReturnMatrix(lines, columns);
            std::vector<std::future<void>> futures;
            for (int i = 0; i < lines; i += blockSize)
            {
                futures.emplace_back(std::async(std::launch::async, [i, blockSize, &ReturnMatrix, this, &Second]()
                    {
                        int end = std::min(i + blockSize, lines);
                        for (int row = i; row < end; ++row)
                        {
                            for (int j = 0; j < columns; ++j)
                            {
                                ReturnMatrix.GetMatrix()[row][j] = matrix[row][j] - Second.GetMatrix()[row][j];
                            }
                        }
                    }));
            }

            for (auto& future : futures)
            {
                future.get();
            }
            return ReturnMatrix;
        }
        else {
            std::cerr << "Error! The matrices must be of the same dimension. Returning a zero matrix with the same dimension as the first one" << std::endl;
            return MatrixSum{ lines, columns };
        }
    }
    template <typename T1>
    auto operator * (Matrix<T1> Second)-> Matrix<decltype(matrix[0][0]*Second.GetMatrix()[0][0])>
    {
        using MatrixMult = Matrix<decltype(matrix[0][0]*Second.GetMatrix()[0][0])>;
        using Mult=decltype(matrix[0][0]*Second.GetMatrix()[0][0]);
        Mult element;
        if (columns == Second.GetLines())
        {
            MatrixMult NewMatrix(lines, Second.GetColumns());
            for (int i = 0; i < lines; i++) {
                for (int j = 0; j < Second.GetColumns(); j++)
                {
                    element = 0;
                    for (int k = 0; k < columns; k++)
                    {
                        element += matrix[i][k] * Second.GetMatrix()[k][j];
                    }
                    NewMatrix.GetMatrix()[i][j] = element;
                }
            }
            return MatrixMult{ lines,Second.GetColumns(),NewMatrix.GetMatrix() };
        }
        else {
            std::cerr << "Error! The matrices must be of the appropriate size. Return a zero matrix with the same dimension as the first one"<<std::endl;
            return MatrixMult{ lines,columns };
        }
    }
    template <typename T1>
    auto Parallel_mult (Matrix<T1> Second)-> Matrix<decltype(matrix[0][0] * Second.GetMatrix()[0][0])>
    {
        using MatrixMult = Matrix<decltype(matrix[0][0] * Second.GetMatrix()[0][0])>;
        using Mult = decltype(matrix[0][0] * Second.GetMatrix()[0][0]);
        Mult element;
        if (columns == Second.GetLines())
        {
            MatrixMult NewMatrix(lines, Second.GetColumns());
            std::vector<std::thread> threads;
            for (int i = 0; i < lines; i++) {
                threads.emplace_back([this, &Second, &NewMatrix, i, &element]() 
                {
                    for (int j = 0; j < Second.GetColumns(); j++)
                    {
                            element = 0;
                            for (int k = 0; k < columns; k++)
                            {
                                element += matrix[i][k] * Second.GetMatrix()[k][j];
                            }
                            NewMatrix.GetMatrix()[i][j] = element;
                    }
                });
            }
            for (std::thread& t : threads) {
                t.join();
            }
            return MatrixMult{ lines,Second.GetColumns(),NewMatrix.GetMatrix() };
        }
        else {
            std::cerr << "Error! The matrices must be of the appropriate size. Return a zero matrix with the same dimension as the first one" << std::endl;
            return MatrixMult{ lines,columns };
        }
    }
    template <typename T1>
    auto Parallel_mult(Matrix<T1> Second, int blockSize) -> Matrix<decltype(matrix[0][0] * Second.GetMatrix()[0][0])>
    {
        using MatrixMult = Matrix<decltype(matrix[0][0] * Second.GetMatrix()[0][0])>;
        using Mult = decltype(matrix[0][0] * Second.GetMatrix()[0][0]);
        Mult element;
        if (columns == Second.GetLines())
        {
            MatrixMult NewMatrix(lines, Second.GetColumns());
            std::vector<std::future<void>> futures;
            for (int i = 0; i < lines; i += blockSize) {
                for (int j = 0; j < Second.GetColumns(); j += blockSize) {
                    futures.push_back(std::async(std::launch::async, [this, &Second, &NewMatrix, blockSize, i, j, &element]()
                        {
                            for (int ii = i; ii < std::min(i + blockSize, lines); ii++)
                            {
                                for (int jj = j; jj < std::min(j + blockSize, Second.GetColumns()); jj++)
                                {
                                    element = 0;
                                    for (int k = 0; k < columns; k++)
                                    {
                                        element += matrix[ii][k] * Second.GetMatrix()[k][jj];
                                    }
                                    NewMatrix.GetMatrix()[ii][jj] = element;
                                }
                            }
                        }));
                }
            }
            for (auto& future : futures) {
                future.get();
            }
            return MatrixMult{ lines, Second.GetColumns(),NewMatrix.GetMatrix() };
        }
        else {
            std::cerr << "Error! The matrices must be of the appropriate size. Return a zero matrix with the same dimension as the first one" << std::endl;
            return MatrixMult{ lines, columns };
        }
    }
    template <typename T1>
    auto operator * (T1 a) -> Matrix<decltype(a*matrix[0][0])>
    {
        using MatrixMult = Matrix<decltype(a*matrix[0][0])>;
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                matrix()[i][j] *= a;
            }
        }
        return MatrixMult{ lines,columns,matrix };
    }
    template <typename T1>
    auto Parallel_mult_numb (T1 a) -> Matrix<decltype(a* matrix[0][0])>
    {
        using MatrixMult = Matrix<decltype(a* matrix[0][0])>;
        std::vector<std::thread> threads; 
        for (int i = 0; i < lines; ++i)
        {
            threads.emplace_back(std::thread([this, a, i]() 
            {
                for (int j = 0; j < columns; ++j)
                {    
                    matrix[i][j] *= a;
                }
            }));
        }

        for (auto& thread : threads)
        {
            thread.join();
        }

        return MatrixMult{ lines, columns, matrix };
    }
    template <typename T1>
    auto Parallel_mult_numb(T1 a, int block_size) -> Matrix<decltype(a* matrix[0][0])>
    {
        using MatrixMult = Matrix<decltype(a* matrix[0][0])>;
        std::vector<std::future<void>> futures;

        for (int i = 0; i < lines; i += block_size)
        {
            futures.push_back(std::async(std::launch::async, [this, a, i, block_size]()
                {
                    for (int x = i; x < std::min(i + block_size, lines); ++x)
                    {
                        for (int j = 0; j < columns; ++j)
                        {
                            matrix[x][j] *= a;
                        }
                    }
                }));
        }

        for (auto& future : futures)
        {
            future.get();
        }

        return MatrixMult{ lines, columns, matrix };
    }
    bool operator != (Matrix<T> Second)
    {
        if (lines == Second.lines && columns == Second.columns)
        {
            for (int i = 0; i < lines; ++i)
            {
                for (int j = 0; j < columns; ++j)
                {
                    if (matrix[i][j] != Second.matrix[i][j])
                    {
                        return true;
                    }
                }
            }
            return false;
        }
        else {
            return true;
        }
    }
    template <typename T1>
    bool operator == (Matrix<T1> Second)
    {
        if (lines == Second.GetLines() && columns == Second.GetColumns())
        {
            for (int i = 0; i < lines; ++i)
            {
                for (int j = 0; j < columns; ++j)
                {
                    if (matrix[i][j] != Second.GetMatrix()[i][j])
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        else {
            return false;
        }
    }
    template<typename T1>
    bool operator ==(T1 a)
    {
        if (columns!=lines) {return false;}
        Matrix<T1> NewMatrix(lines, columns);
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                if (i == j) {
                    NewMatrix.GetMatrix()[i][i] = a;
                }
            }
        }
        if (*this==NewMatrix)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    template<typename T1>
    bool operator !=(T1 a)
    {
        if (columns!=lines) {return false;}
        Matrix<T1> NewMatrix(lines, columns);
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                if (i == j) {
                    NewMatrix.GetMatrix()[i][i] = a;
                }
            }
        }
        if (*this!=NewMatrix)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    Matrix<double> operator ! ()
    {
        Matrix<double> AlgebraicAdditions(lines, columns);
        T det = this->Determinant();
        if (det == 0||(lines!=columns))
        {
            throw "Error! The matrix must be square and have a non-zero determinant.";
            return { lines,columns };
        }
        for (int i = 0; i < lines; ++i) {
            for (int j = 0; j < lines; ++j) {
                Matrix Addition(*this, i, j);
                AlgebraicAdditions.GetMatrix()[i][j] = pow(-1, i + j) * Addition.Determinant() / det;
            }
        }
        AlgebraicAdditions.Transposition();
        return Matrix<double>{ AlgebraicAdditions.GetLines(),AlgebraicAdditions.GetColumns(),AlgebraicAdditions.GetMatrix() };
    }
    Matrix operator * (int a)
    {
        Matrix ReturnMatrix(lines, columns);
        for (int i = 0; i < lines; ++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                ReturnMatrix.matrix[i][j]=matrix[i][j] * a;
            }
        }
        return ReturnMatrix;
    }
    Matrix operator =(Matrix Second)
    {
        lines=Second.GetLines();
        columns=Second.GetColumns();
        for (int i = 0; i < lines; ++i)
        {
            delete[] matrix[i];
        }
        delete[] matrix;
        T** newMatrix=new T*[Second.GetLines()];
        for (int i=0;i<Second.GetLines();++i)
        {
            newMatrix[i]=new T[columns];
        }
        for (int i=0;i<lines;++i)
        {
            for (int j = 0; j < columns; ++j)
            {
                matrix[i][j] = Second.matrix[i][j];
            }
        }
        this->matrix = newMatrix;
    }
};

