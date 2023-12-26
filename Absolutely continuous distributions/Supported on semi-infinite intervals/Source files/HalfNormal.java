public class HalfNormal extends MathLibrary{
    double sigma;

    HalfNormal(double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.sigma = sigma;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (root(2, 2) / (sigma * root(PI, 2))) * euler(-1 * power(x, 2) / (2 * power(sigma, 2)));
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return errorFunction(x / (sigma * root(2, 2)));
    }

    double Quantile(double x) {
        return (sigma / root(2, 2)) * INVERF(x);
    }

    double Mean() {
        return (sigma * root(2, 2)) / root(PI, 2);
    }

    double Median() {
        return sigma * root(2, 2) * INVERF(0.5);
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        return power(sigma, 2) * (1 - (2 / PI));
    }

    double Skewness() {
        return (root(2, 2) * (4 - PI)) / ((PI - 2) * root(PI - 2, 2));
    }

    double EXkurtosis() {
        return (8 * (PI - 3)) / power(PI - 2, 2);
    }

    double Entropy() {
        return 0.5 * log(2 * PI * euler(1) * power(sigma, 2), 2) - 1;
    }
}
