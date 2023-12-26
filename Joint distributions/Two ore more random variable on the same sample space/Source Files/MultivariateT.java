public class MultivariateT extends MathLibrary{
    double[][] mu;
    double[][] sigma;
    double nu;

    MultivariateT(double[][] mu, double[][] sigma, double nu) {
        if (mu.length > 1) throw new IllegalArgumentException("mu = double[0][x]");
        if (sigma.length != sigma[0].length) throw new IllegalArgumentException("sigma = double[x][x]");
        if (mu[0].length != sigma.length) throw new IllegalArgumentException("mu[0].length = sigma.length");
        if (nu <= 0) throw new IllegalArgumentException("nu > 0");
        this.mu = mu;
        this.sigma = sigma;
        this.nu = nu;
    }

    double[][] PDF(double[][] x) {
        double part1 = gamma((nu + mu.length) / 2);
        double part2 = gamma(nu / 2) * doublePower(nu, mu.length / 2) * doublePower(PI, mu.length / 2) * root(determinant(sigma), 2);
        double[][] part3 = matrixMultiplication(matrixByMatrix(transpose(matrixSubtraction(x, mu)), transpose(matrixByMatrix(inverseMatrix(sigma), transpose(matrixSubtraction(x, mu))))), 1 / nu);
        return matrixMultiplication(matrixPower(part3, (int) (-1 * (nu + mu.length) / 2)), part1 / part2);
    }

    double[][] Mean() {
        if (nu > 1) return mu;
        else return new double[0][0];
    }

    double[][] Median() {
        return mu;
    }

    double[][] Mode() {
        return mu;
    }

    double[][] Variance() {
        if (nu > 2) {
            return matrixMultiplication(sigma, nu / (nu - 2));
        } else {
            return new double[0][0];
        }
    }

    int Skewness() {
        return 0;
    }
}
