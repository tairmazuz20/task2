# Copy matrix#
def copyMatrix(mat, size):
    newMat = []
    for i in range(size):
        rowList = []
        for j in range(size):
            rowList.append(mat[i][j])
        newMat.append(rowList)
    return newMat


# Check if the diagonal is dominant#
def isDDM(matrix, n):
    mat = copyMatrix(matrix, n)
    for i in range(0, n):
        sum1 = 0
        for j in range(0, n):
            sum1 = sum1 + abs(mat[i][j])
        sum1 = sum1 - abs(mat[i][i])
        if abs(mat[i][i]) < sum1:
            return False
    return True


# Matrix construction D#
def initD(matrix, n):
    mat = copyMatrix(matrix, n)
    for i in range(n):
        num = abs(mat[i][i])
        for j in range(n):
            (mat[i][j]) = 0
        (mat[i][i]) = num
    return mat


# Matrix multiplication 3X3#
def multiMatrix(m1, m2):
    res = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(len(m1)):
        # iterate through columns of Y
        for j in range(len(m2[0])):
            # iterate through rows of Y
            for k in range(len(m2)):
                res[i][j] += m1[i][k] * m2[k][j]
    return res


# Matrix multiplication 3X1#
def multiM(m1, m2):
    res1 = [[0, ], [0, ], [0, ]]
    for i in range(len(m1)):
        # iterate through columns of Y
        for j in range(len(m2[0])):
            # iterate through rows of Y
            for k in range(len(m1)):
                res1[i][j] += (m1[i][k] * m2[k][j])

    return res1


def getMatrixMinor(m, i, j):
    return [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]


# Determination of determination#
def getMatrixDet(m):
    # base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1) ** c) * m[0][c] * getMatrixDet(getMatrixMinor(m, 0, c))
    return determinant


# Arrangement of a matrix#
def orderMatrix(a):
    for i in range(len(a)):
        for j in range(len(a)):
            if i == j and a[i][j] == 0:
                q = a[i]
                for n in range(len(a)):
                    if a[n][i] != 0:
                        a[i] = a[n]
                        a[n] = q


# Inversion matrix calculation#
def reverseMatrix(m1):
    if getMatrixDet(m1) != 0:
        orderMatrix(m1)
        IA = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        res = multiMatrix(IA, IA)
        for j in range(len(m1)):
            for i in range(len(m1)):
                if i >= j:
                    IA = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                    if i == j:
                        dev = m1[i][j]
                        for k in range(len(m1)):
                            m1[i][k] /= dev
                        IA[i][j] /= dev
                        res = multiMatrix(IA, res)
                    else:
                        if m1[i][j] != 0:
                            mu = -m1[i][j]
                            for k in range(len(m1)):
                                m1[i][k] += mu * m1[j][k]
                            IA[i][j] = mu
                            res = multiMatrix(IA, res)
        for j in range(len(m1) - 1, -1, -1):
            for i in range(len(m1) - 2, -1, -1):
                if j > i:
                    IA = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                    if m1[i][j] != 0:
                        mu = -m1[i][j]
                        for k in range(len(m1)):
                            m1[i][k] += mu * m1[j][k]
                        IA[i][j] = mu
                        res = multiMatrix(IA, res)
        return res


# Finding Matrix L and U#
def LU(mat, n):
    ll = [[0, 0, 0],
          [0, 0, 0],
          [0, 0, 0]]
    U = [[0, 0, 0],
         [0, 0, 0],
         [0, 0, 0]]
    for i in range(len(mat)):
        for j in range(len(mat)):
            if i != j:
                if i > j:
                    ll[i][j] = mat[i][j]
                else:
                    U[i][j] = mat[i][j]

    return ll, U


# The minus matrix#
def minusMat(mat):
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            mat[i][j] = -mat[i][j]
    return mat


# Connecting matrices#
def connectingMat(a, b):
    for i in range(len(a)):
        for j in range(len(a[0])):
            a[i][j] += b[i][j]
    return a


# Finding matrices L and U and D and the inverse D#
def solve(mat, n):
    L, U = LU(mat, n)
    D = initD(mat, n)
    cD = copyMatrix(D, n)
    D1 = reverseMatrix(cD)
    return L, U, D, D1


# Subtraction of matrices#
def subMat(m1, m2):
    subM = 0
    for i in range(len(m1)):
        for j in range(1):
            subM += m1[i][j] - m2[i][j]
    return subM


# Solution using the Yak-obi method#
def yakovi(mat, n, B):
    print("yakovi")
    Xr = [[0], [0], [0]]
    Xr1 = [[0], [0], [0]]
    L, U, D, D1 = solve(mat, n)
    H = D1
    G = minusMat(multiMatrix(D1, connectingMat(L, U)))
    Xr1 = connectingMat(multiM(G, Xr), multiM(H, B))
    while abs(subMat(Xr1, Xr)) > 0.001:
        Xr = Xr1
        Xr1 = connectingMat(multiM(G, Xr), multiM(H, B))
        print(Xr1)
    return Xr1


# Solution using the Seizel method#
def Seizel(mat, n, B):
    print("Seizel")
    Xr = [[0], [0], [0]]
    Xr1 = [[0], [0], [0]]
    L, U, D, D1 = solve(mat, n)
    H = reverseMatrix(connectingMat(L, D))
    G = minusMat(multiMatrix(H, U))
    Xr1 = connectingMat(multiM(G, Xr), multiM(H, B))
    print(Xr1)
    while abs(subMat(Xr1, Xr)) > 0.001:
        Xr = Xr1
        Xr1 = connectingMat(multiM(G, Xr), multiM(H, B))
        print(Xr1)
    return Xr1


def main():
    n = 3
    mat = [[4, 2, 0],
           [2, 10, 4],
           [0, 4, 5]]
    B = [[2],
         [6],
         [5]]
    dd = initD(mat, n)
    if not isDDM(mat, n):
        print("In the matrix there is no dominant diagonal")
    else:
        yakovi(mat, n, B)
        Seizel(mat, n, B)


main()
