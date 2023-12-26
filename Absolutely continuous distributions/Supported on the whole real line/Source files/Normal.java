public class Normal extends MathLibrary{
    double mu;
    double sigma;

    Normal(double mu, double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
    }

    double PDF(double x) {
        double exponent = -0.5 * power((x - mu) / root(sigma, 2), 2);
        return (1 / root(2 * PI * sigma, 2)) * euler(exponent);
    }

    double CDF(double x) {
        return 0.5 * (1 + errorFunction((x - mu) / root(2 * sigma, 2)));
    }

    double Quantile(double p) {
        return mu + (root(2 * sigma, 2) * INVERF(2 * p - 1));
    }

    double Mean() {
        return mu;
    }

    double Median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return sigma;
    }

    double MAD() {
        return root(2 * sigma / PI, 2);
    }

    int Skewness() {
        return 0;
    }

    int EXkurtosis() {
        return 0;
    }

    double Entropy() {
        return 0.5 * ln(2 * PI * euler(1) * sigma);
    }

    double MGF(double t) {
        return euler(mu * t + (sigma * power(t, 2) / 2));
    }

    double[][] FisherInformation() {
        double[][] values = {
            {1 / sigma, 0},
            {0, 1 / (2 * power(sigma, 2))}
        };
        return values;
    }

    double KullbackLeiblerDivergence(double mu2, double sigma2) {
        double part1 = power(root(sigma2, 2) / root(sigma2, 2), 2) + (power(mu - mu2, 2) / sigma);
        return 0.5 * (part1 - 1 + ln(sigma2 / sigma));
    }
}
