public class NegativeMultinomial extends MathLibrary{
    int x0;
    double[][] p;
    double p0;

    NegativeMultinomial(int x0, double[][] p) {
        if (x0 <= 0) throw new IllegalArgumentException("x0 > 0");
        if (p.length > 1) throw new IllegalArgumentException("p = double[0][x]");
        this.x0 = x0;
        this.p = p;
        double value = 0;
        for (int i = 0; i < p.length; i++) {
            value += p[0][i];
        }
        this.p0 = 1 - value;
    }

    double PMF(int[] x) {
        if (x.length < 1 || x.length > p.length) throw new IllegalArgumentException("1 <= x.length <= p.length");
        double value1 = 0;
        double value2 = 1;
        for (int i = 0; i < x.length; i++) {
            value1 += x[i];
        }
        for (int i = 0; i < x.length; i++) {
            value2 *= power(p[0][i], x[i]) / factorial(x[i]);
        }
        return gamma(value1) * ((power(p0, x0) / gamma(x0)) * value2);
    }

    double[][] Mean() {
        return matrixMultiplication(p, x0 / p0);
    }

    double[][] Variance() {
        double[][] part1 = matrixMultiplication(matrixByMatrix(p, transpose(p)), x0 / power(p0, 2));
        double[][] part2 = matrixMultiplication(digaonal(p), x0 / p0);
        return matrixSubtraction(part1, matrixMultiplication(part2, -1));
    }

    double MGF(double[] t) {
        double value = 0;
        for (int i = 0; i < p.length; i++) {
            value += p[0][i] * euler(t[i]);
        }
        return power(p0 / (1 - value), x0);
    }
}
