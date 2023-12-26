public class MathLibrary {
    protected double PI = 3.1415;
    protected double EulerMascheroni = 0.5772;

    protected int round(double x) {
        if (x - (int) x < 0.5) {
            return (int) x;
        } else {
            return (int) x + 1;
        }
    }

    protected double absolute(double value) { // zwraca wartość bezwzgęldną z liczby
        if (value >= 0) {
            return value;
        } else {
            return value * -1;
        }
    }

    protected double max(double x1, double x2) {
        if (x1 > x2) return x1;
        else return x2;
    }

    protected double min(double x1, double x2) {
        if (x1 < x2) return x1;
        else return x2;
    }

    protected int sign(double x) {
        if (x < 0) return -1;
        else if (x == 0) return 0;
        else return 1;
    }
    
    protected double power(double base, int exponent) {
        double value = 1;
        if (exponent < 0) {
            base = 1 / base;
            exponent *= -1;
        }
        if (base != 0) {
            for (int i = 0; i < exponent; i++) {
                value *= base;
            }
        }
        return value;
    }

    protected double doublePower(double base, double exponent) {
        int integerPart = (int) exponent;
        double fractionalPart = exponent - integerPart;
        double integerResult = power(base, integerPart);
        double fractionalResult = 1;
        if (fractionalPart != 0) {
            fractionalResult = euler(fractionalPart * ln(base));
        }
        return integerResult * fractionalResult;
    }

    protected double ln(double x) {
        if (x < 1) return -ln(1 / x);
        if (x == 1) return 0;

        double term = (x - 1) / (x + 1);
        double term2 = term * term;
        double result = term;
        int i = 1;
        while (absolute(term) > 1e-10) { // precision of 10 decimal places
            term *= term2;
            i += 2;
            result += term / i;
        }
        return 2 * result;
    }

    protected double log(double x, int podstawa) {
        return ln(x) / ln(podstawa);
    }

    protected double polylogarithm(double x, double s) {
        double value = 0;
        for (int i = 1; i < 9; i++) {
            value += power(x, i) / doublePower(i, s);
        }
        return value;
    }

    protected double root(double x, int n) {
        if (n == 0)
            return 0;
        double value = 2;
        for (int i = 0; i < 15; i++) {
            value = (((double) n - 1) / n * value) + (((double) x / n) * (1.0 / power(value, n - 1)));
        }
        return value;
    }

    protected long factorial(int x) {
        long value = 1;
        if (x != 0) {
            for (int i = x; i > 0; i--) {
                value *= i;
            }
        }
        return value;
    }

    protected double cosine(double x) {
        double value = 0;
        for (int n = 0; n < 11; n++) {
            value += power(-1, n) * (power(x, 2 * n) / (double) factorial(2 * n));
        }
        return value;
    }

    protected double sine(double x) {
        double value = 0;
        for (int n = 0; n < 11; n++) {
            value += power(-1, n) / factorial(2 * n + 1) * power(x, 2 * n + 1);
        }
        return value;
    }

    protected double tangent(double x) {
        return sine(x) / cosine(x);
    }

    protected double hyperbolicSine(double x) {
        return (euler(x) - euler(-1 * x)) / 2;
    }

    protected double hyperbolicCosine(double x) {
        double value = 0;
        for (int n = 0; n < 7; n++) {
            value += (power(x, 2 * n)) / factorial(2 * n);
        }
        return value;
    }

    protected double hyperbolicTangent(double x) {
        return (double) (euler(2 * x) - 1) / (euler(2 * x) + 1);
    }

    protected double arcSine(double x) {
        double value = 0;
        for (int n = 0; n < 7; n++) {
            double part = factorial(2 * n)  / (power(4, n) * power(factorial(n), 2) * (2 * n + 1));
            value += part * power(x, 2 * n + 1);
        }
        return value;
    }

    protected double arcTangent(double x) {
        double value = 0;
        for (int n = 0; n < 7; n++) {
            value += power(-1, n) / (2 * n + 1) * power(x, 2 * n + 1);
        }
        return value;
    }

    protected double hyperbolicArcSine(double x) {
        return ln(x + root(power(x, 2) + 1, 2));
    }

    protected double hyperbolicArcTangent(double x) {
        return 0.5 * ln((1 + x) / (1 - x));
    }

    protected double euler(double x) {
        double value = 1;
        if (x != 0) {
            for (int n = 50; n > 0; --n) {
                value = 1 + x * value / n;
            }
        }
        return value;
    }

    protected double binomial(double n, double k) {
        return gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1));
    }

    protected double pochhammer(double x, int n) {
        double value = 1;
        for (int i = 0; i < n - 1; i++) {
            if (x + i == 0) value *= 1;
            else value *= x + i;
        }
        return value;
    }

    protected double gamma(double x) {
        if (x <= 2) return 1.0;
        double temp = (x - 0.5) * ln(x + 4.5) - (x + 4.5);
        double CD = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
        + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
        +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
        double logGamma = temp + ln(CD * root(2 * PI, 2));
        return euler(logGamma);
    }

    protected double incompleteLowerGamma(double s, double x) {
        double value = 0;
        for (int i = 0; i < s; i++) {
            value += power(x, i) / factorial(i);
        }
        return gamma(s) * (1 - (euler(-1 * x) * value));
    }

    protected double incompleteUpperGamma(double s, double x) {
        return gamma(s) - incompleteLowerGamma(s, x);
    }

    protected double regularizedGamma(double s, double x) {
        return incompleteLowerGamma(s, x) / gamma(s);
    }

    protected double polygamma(int m, double x) {
        double value = 0;
        for (int i = 0; i < 7; i++) {
            value += 1 / power(x + i, m + 1);
        }
        return power(-1, m + 1) * factorial(m) * value;
    }

    protected double digamma(double z) {
        return ln(z) - (1 / (2 * z));
    }

    protected double trigamma(double x) {
        double value = 0;
        for (int i = 0; i < 30; i++) {
            value += 1 / power(x + i, 2);
        }
        return value;
    }

    protected double beta(double x, double y) {
        return (gamma(x) * gamma(y)) / gamma(x + y);
    }

    protected double incompleteBeta(double a, double b, double z) {
        double value = 0;
        for (int i = 0; i < 10; i++) {
            System.out.println(pochhammer(1 - b, i));
            value += pochhammer(1 - b, i) / (factorial(i) * (a + i)) * power(z, i);
        }
        return doublePower(z, a) * value;
    }

    protected double regularizedIncompleteBeta(double a, double b, double z) {
        return incompleteBeta(a, b, z) / beta(a, b);
    }

    protected double errorFunction(double x) {
        return hyperbolicTangent((x * PI) / root(6, 2));
    }

    protected double INVERF(double x) {
        double z;
        double a  = 0.147;                                                   
        int the_sign_of_x = sign(x);

        if(0 != x) {
            double ln_1minus_x_sqrd = ln(1-x*x);
            double ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
            double ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
            double ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2/(PI * a));
            double first_sqrt = root((ln_etc_by2_plus2*ln_etc_by2_plus2)-ln_1minusxx_by_a, 2);
            double second_sqrt = root(first_sqrt - ln_etc_by2_plus2, 2);
            z = second_sqrt * the_sign_of_x;
        } else { // x is zero
            z = 0;
        }
  return z;
    }

    protected double ERFC(double x) {
        return 1 - errorFunction(x);
    }

    private double lambertW(double x) {
        double w = ln(x);
        for (int i = 0; i < 10; i++) {
            double f = w * euler(w) - x;
            double df = (w + 1) * euler(w);
            w -= f / df;
        }
        return w;
    }

    protected double INVERFC(double x) {
        return root(lambertW(2 / (PI * power(x, 2))), 2) / root(2, 2) - (0.5 * x);
    }

    protected double zeta(double s) {
        double value = 0;
        for (int i = 1; i < 20; i++) {
            value += 1.0 / doublePower(i, s);
        }
        return value;
    }

    protected double generalizedHarmonicNumber(int n, double m) {
        double value = 0;
        for (int i = 1; i <= n; i++) {
            value += 1.0 / doublePower(i, m);
        }
        return value;
    }

    protected double logit(double x) {
        return ln(x / (1 - x));
    }

    protected double modifiedBessel(double alpha, double x) {
        double value = 0;
        for (int m = 0; m < 20; m++) {
            value += (1 / (factorial(m) * gamma(m + alpha + 1))) * doublePower(x / 2, 2 * m + alpha);
        }
        return value;
    }

    protected double modifiedBessel2(double nu, double z) {
        double a = 0.0;
        double b = ln(2 * PI * absolute(z) / hyperbolicTangent(PI * nu));
        double stepSize = (b - a) / 10000;

        double integral = 0.0;
        for (double x = a; x <= b; x += stepSize) {
            integral += euler(-z * hyperbolicCosine(x)) * hyperbolicCosine(nu * x) * stepSize;
        }

        return integral;
    }

    protected double eulerFunction(double q) {
        double value = 1;
        for (int i = 1; i < 15; i++) {
            value *= (1 - power(q, i));
        }
        return value;
    }

    protected double marcumQ(double M, double a, double b) {
        double value = 0;
        for (int k = 0; k < 7; k++) {
            value += (1 / factorial(k)) * (incompleteLowerGamma(M + k, power(b, 2) / 2) / gamma(M + k)) * power(power(a, 2) / 2, k);
        }
        return 1 - (euler(-1 * power(a, 2) / 2) * value);
    }

    protected double laguerre(double bottom, double up, double x) {
        double value = 0;
        for (int i = 0; i < (int) (bottom + 1); i++) {
            value += power(-1, i) * binomial(bottom + up, bottom - i) * (power(x, i) / factorial(i));
        }
        return value;
    }

    protected double hypergeometric(double a, double b, double c, double z) {
        double value = 0;
        for (int i = 0; i < 20; i++) {
            value += ((pochhammer(a, i) * pochhammer(b, i)) / pochhammer(c, i)) * (power(z, i) / factorial(i));
        }
        return value;
    }

    protected double exponentialIntergral(double x) {
        double xn = -x;
        double Sn = -x;
        double Sm1 = 0.0;
        double hsum = 1.0;
        double g = 0.5772156649015328606065121;
        double y = 1.0;
        double factorial = 1.0;
 
        while (absolute(Sn - Sm1) > 10.0 * Math.ulp(1.0) * absolute(Sm1)) {
            Sm1 = Sn;
            y += 1.0;
            xn *= (-x);
            factorial *= y;
            hsum += (1.0 / y);
            Sn += hsum * xn / factorial;
        }
        return (g + ln(absolute(x)) - euler(x) * Sn);
    }

    protected double logarithmicIntegral(double x) {
        return exponentialIntergral(ln(x));
    }

    protected double[][] matrixSubtraction(double[][] x, double[][] y) {
        double[][] values = new double[x.length][x[0].length];
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < x[0].length; j++) {
                values[i][j] = x[i][j] - y[i][j];
            }
        }
        return values;
    }

    protected double[][] matrixMultiplication(double[][] x, double value) {
        double[][] values = new double[x.length][x[0].length];
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < x[0].length; j++) {
                values[i][j] = x[i][j] * value;
            }
        }
        return values;
    }

    protected double[][] matrixByMatrix(double[][] a, double[][] b) {
        int aRows = a.length;
        int aColumns = a[0].length;
        int bRows = b.length;
        int bColumns = b[0].length;

        if (aColumns != bRows) {
            throw new IllegalArgumentException("Number of columns of first matrix must be equal to number of rows of second matrix.");
        }
 
        double[][] result = new double[aRows][bColumns];
 
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                for (int k = 0; k < aColumns; k++) {
                   result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
 
        return result;
    }

    protected double[][] matrixPower(double[][] a, int exponent) {
        int n = a.length;
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = (i == j ? 1 : 0);
            }
        }
        while (exponent > 0) {
            if ((exponent & 1) != 0) {
                result = matrixByMatrix(result, a);
            }
            a = matrixByMatrix(a, a);
            exponent >>= 1;
        }
  
        return result;
    }

    protected double[][] eulerMatrix(double[][] x) {
        double[][] values = new double[x.length][x.length];

        for (int i = 0; i < 20; i++) {
            values = matrixSubtraction(values, matrixMultiplication(matrixMultiplication(matrixPower(x, i), 1.0 / factorial(i)), -1));
        }
        return values;
    }

    protected double[][] transpose(double[][] martrix) {
        double[][] values = new double[martrix[0].length][martrix.length];
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                values[i][j] = martrix[j][i];
            }
        }
        return values;
    }

    private double[][] generateSubArray(double A[][] , int j1){
        double[][] m = new double[A.length - 1][];
        for (int k = 0; k < (A.length - 1); k++) {
            m[k] = new double[A.length - 1];
        }
        for (int i = 1; i < A.length; i++){
              int j2=0;
              for (int j = 0; j < A.length; j++){
                  if(j == j1)
                        continue;
                  m[i - 1][j2] = A[i][j];
                  j2++;
              }
        }
        return m;
    }

    protected double determinant(double A[][]){
        double res;
        double[][] m;

        if (A.length == 1) res = A[0][0];
        else if (A.length == 2) res = A[0][0]*A[1][1] - A[1][0]*A[0][1];
        else{
            res=0;
            for (int j1 = 0; j1 < A.length; j1++){
                m = generateSubArray(A, j1);
                res += Math.pow(-1.0, 1.0 + j1 + 1.0) * A[0][j1] * determinant(m);
            }
        }
        return res;
    }

    private double[][] submatrix(double[][] matrix, int row, int column) {
        double[][] submatrix = new double[matrix.length - 1][matrix.length - 1];

        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; i != row && j < matrix[i].length; j++)
                if (j != column)
                    submatrix[i < row ? i : i - 1][j < column ? j : j - 1] = matrix[i][j];
        return submatrix;
    }

    protected double[][] inverseMatrix(double[][] matrix) {
        double[][] inverse = new double[matrix.length][matrix.length];

        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                inverse[i][j] = Math.pow(-1, i + j)
                        * determinant(submatrix(matrix, i, j));

        double det = 1.0 / determinant(matrix);
        for (int i = 0; i < inverse.length; i++) {
            for (int j = 0; j <= i; j++) {
                double temp = inverse[i][j];
                inverse[i][j] = inverse[j][i] * det;
                inverse[j][i] = temp * det;
            }
        }

        return inverse;
  }

    protected double trace(double[][] x) {
        double value = 0;
        for (int i = 0; i < x.length; i++) {
            value += x[i][i];
        }
        return value;
    }

    protected double[][] digaonal(double[][] x) {
        double[][] result = new double[x.length][x[0].length];
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < x[i].length; j++) {
                if (i == j) result[i][j] = x[i][j];
                else result[i][j] = 0;
            }
        }
        return result;
    }
}