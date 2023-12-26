public class LogNormal extends MathLibrary{
    double mu;
    double sigma;

    LogNormal(double mu, double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = power(sigma, 2);
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part = euler(-1 * power(ln(x) - mu, 2) / (2 * sigma));
        return (1 / (x * root(2 * PI * sigma, 2))) * part;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return 0.5 * (1 + errorFunction((ln(x) - mu) / root(2 * sigma, 2)));
    }

    double Quantile(double p) {
        return euler(mu + (root(2 * sigma, 2) * INVERF(2 * p - 1)));
    }

    double Mean() {
        return euler(mu + (sigma / 2));
    }

    double Median() {
        return euler(mu);
    }

    double Mode() {
        return euler(mu - sigma);
    }

    double Variance() {
        return (euler(sigma) - 1) * euler(2 * mu + sigma);
    }

    double Skewness() {
        return (euler(sigma) + 2) * root(euler(sigma) - 1, 2);
    }

    double EXkurtosis() {
        return euler(4 * sigma) + (2 * euler(3 * sigma)) + (3 * euler(2 * sigma)) - 6;
    }

    double Entropy() {
        return log(root(2 * PI * sigma, 2) * euler(mu + 0.5), 2);
    }

    double[][] FisherInformation() {
        double[][] values = {
            {1 / sigma, 0},
            {0, 2 / sigma}
        };
        return values;
    }

    double ExpectedShortfall(double p) {
        double part1 = ERFC((sigma / root(2, 2)) - INVERF(2 * p - 1)) / 2;
        double part2 = (1 - p) * euler(mu + (sigma / 2));
        return part1 * part2;
    }
}
