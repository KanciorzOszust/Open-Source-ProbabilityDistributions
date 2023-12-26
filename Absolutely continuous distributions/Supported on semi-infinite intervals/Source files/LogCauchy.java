public class LogCauchy extends MathLibrary{
    double mu;
    double sigma;

    LogCauchy(double mu, double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part = sigma / (power(ln(x) - mu, 2) + power(sigma, 2));
        return (1 / (x * PI)) * part;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return (1 / PI) * arcTangent((ln(x) - mu) / sigma) + 0.5;
    }

    double Mean() {
        return Double.POSITIVE_INFINITY;
    }

    double Median() {
        return euler(mu);
    }

    double Variance() {
        return Mean();
    }
}
