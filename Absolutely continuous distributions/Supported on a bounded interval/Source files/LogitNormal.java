public class LogitNormal extends MathLibrary{
    double sigma;
    double mu;

    LogitNormal(double sigma, double mu) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.sigma = sigma;
        this.mu = mu;
    }

    double PDF(double x) {
        if (x <= 0 || x >= 1) throw new IllegalArgumentException("0 < x < 1");
        double exponent = -1 * power(logit(x) - mu, 2) / (2 * sigma);
        return 1.0 / (sigma * root(2 * PI, 2)) * euler(exponent) * (1 / (x * (1 - x)));
    }

    double CDF(double x) {
        if (x <= 0 || x >= 1) throw new IllegalArgumentException("0 < x < 1");
        double argument = (logit(x) - mu) / root(2 * sigma, 2);
        return 0.5 * (1 + errorFunction(argument));
    }
}
