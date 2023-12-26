public class MultivariateNormal extends MathLibrary{
    double[][] mu;
    double[][] sigma;

    MultivariateNormal(double[][] mu, double[][] sigma) {
        if (mu.length > 1) throw new IllegalArgumentException("mu = double[0][x]");
        if (sigma.length != sigma[0].length) throw new IllegalArgumentException("sigma = double[x][x]");
        if (mu[0].length != sigma.length) throw new IllegalArgumentException("mu[0].length = sigma.length");
        this.mu = mu;
        this.sigma = sigma;
    }

    double[][] PDF(double[][] x) {
        double[][] exponent = matrixMultiplication(transpose(matrixSubtraction(x, mu)), -0.5);
        exponent = matrixByMatrix(transpose(exponent), matrixByMatrix(inverseMatrix(sigma), transpose(matrixSubtraction(x, mu))));
        return matrixMultiplication(eulerMatrix(exponent), doublePower(2 * PI, -1 * sigma.length / 2.0) * doublePower(determinant(sigma), -1.0 / 2));
    }

    double[][] Mean() {
        return mu;
    }

    double[][] Mode() {
        return mu;
    }

    double[][] Variance() {
        return sigma;
    }

    double Entropy() {
        return (sigma.length / 2) * ln(2 * PI * euler(1)) + (0.5 * ln(determinant(sigma)));
    }

    double[][] MGF(double[][] t) {
        double[][] exponent = matrixByMatrix(transpose(mu), t); 
        double[][] exponent2 = matrixMultiplication(matrixByMatrix(transpose(t), transpose(matrixByMatrix(sigma, transpose(t)))), 0.5);  
        return eulerMatrix(matrixSubtraction(exponent, matrixMultiplication(exponent2, -1)));
    }

    double KullbackLeiblerDivergence(double[][] sigma1) {
        double part1 = trace(matrixByMatrix(inverseMatrix(sigma1), sigma)) - sigma.length;
        double part2 = ln(determinant(sigma1) / determinant(sigma));
        return 0.5 * (part1 + part2);
    }
}
