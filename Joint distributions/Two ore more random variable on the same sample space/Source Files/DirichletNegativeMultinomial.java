public class DirichletNegativeMultinomial extends MathLibrary{
    double x0;
    double alpha0;
    double[][] alpha;
    double pointAlpha;

    DirichletNegativeMultinomial(double x0, double alpha0, double[][] alpha) {
        if (x0 <= 0 || alpha0 <= 0) throw new IllegalArgumentException("x0, alpha0 > 0");
        if (alpha.length > 1) throw new IllegalArgumentException("alpha = double[0][x]");
        for (int i = 0; i < alpha.length; i++) {
            if (alpha[0][i] <= 0) throw new IllegalArgumentException("alpha[0][i] > 0");
            this.pointAlpha += alpha[0][i];
        }
        this.x0 = x0;
        this.alpha = alpha;
        this.alpha0 = alpha0;
    }

    private double pointX(int[] x) {
        double value = 0;
        for (int i = 0; i < x.length; i++) {
            value += x[i];
        }
        return value;
    }

    double PMF(int[] x) {
        if (x.length != alpha.length) throw new IllegalArgumentException("x.length = alpha.length");
        double value = 1;
        for (int i = 0; i < x.length; i++) {
            value *= gamma(x[i] + alpha[0][i]) / (factorial(x[i]) * gamma(alpha[0][i]));
        }
        return (beta(pointX(x), pointAlpha) / beta(x0, alpha0)) * value;
    }

    double[][] Mean() {
        if (alpha0 > 1) {
            System.out.println(x0);
            return matrixMultiplication(alpha, x0 / (alpha0 - 1));
        } else {
            return new double[0][0];
        }
    }

    double[][] Variance() {
        if (alpha0 > 2) {
            double part1 = x0 * (x0 + alpha0 - 1);
            double part2 = power(alpha0 - 1, 2) * (alpha0 - 2);
            double[][] part3 = matrixByMatrix(alpha, transpose(alpha));
            double[][] part4 = matrixMultiplication(digaonal(alpha), alpha0 - 1);
            return matrixMultiplication(matrixSubtraction(part3, matrixMultiplication(part4, -1)), part1 / part2);
        } else {
            return new double[0][0];
        }
    }
}
