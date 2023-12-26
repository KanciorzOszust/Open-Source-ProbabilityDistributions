public class ExponentiallyModifiedGaussian extends MathLibrary{
    double mu;
    double sigma;
    double lambda;

    ExponentiallyModifiedGaussian(double mu, double sigma, double lambda) {
        if (sigma <= 0 || lambda <= 0) throw new IllegalArgumentException("sigma, lambda > 0");
        this.mu = mu;
        this.sigma = sigma;
        this.lambda = lambda;
    }

    double PDF(double x) {
        double part1 = (lambda / 2) * euler((lambda / 2) * (2 * mu + (lambda * sigma) - (2 * x)));
        double part2 = ERFC((mu + (lambda * sigma) - x) / root(2 * sigma, 2));
        return part1 * part2;
    }

    private double PHI(double x) {
        return 0.5 * (1 + errorFunction((x - mu) / root(2 * sigma, 2)));
    }

    double CDF(double x) {
        return PHI(x) - (0.5 * PDF(x));
    }

    double Mean() {
        return mu + (1 / lambda);
    }

    double Mode() {
        double r = 1 / lambda;
        double xm = mu - (root(2 * sigma, 2) * INVERFC((r / root(sigma, 2)) * root(2 / PI, 2))) + (sigma / r);
        double part = -0.5 * power((mu - xm) / root(sigma, 2), 2);
        double h = (1 / root(2 * PI * sigma, 2)) * euler(-0.5 * power(-1 * mu / root(sigma, 2), 2));
        return h * euler(part);
    }

    double Variance() {
        return sigma + (1 / power(lambda, 2));
    }

    double Skewness() {
        return (2 / sigma * root(sigma, 2) * power(lambda, 3)) * doublePower(1 + (1 / (sigma * power(lambda, 2))), -1.5);
    }

    double EXkurtosis() {
        double part1 = 3 * (1 + (2 / (sigma * power(lambda, 2))) + (3 / (power(sigma, 2) * power(lambda, 4))));
        double part2 = power(1 + (1 / (sigma * power(lambda, 2))), 2);
        return part1 / part2 - 3;
    }

    double MGF(double t) {
        double part = euler(mu * t + (0.5 * sigma * power(t, 2)));
        return (1 / (1 - (t / lambda))) * part;
    }
}
