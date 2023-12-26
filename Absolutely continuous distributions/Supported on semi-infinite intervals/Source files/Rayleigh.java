public class Rayleigh extends MathLibrary{
    double sigma;

    Rayleigh(double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.sigma = power(sigma, 2);
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (x / sigma) * euler(-1 * power(x, 2) / (2 * sigma));
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - euler(-1 * power(x, 2) / (2 * sigma));
    }

    double Quantile(double F) {
        return root(sigma * -2 * ln(1 - F), 2);
    }

    double Mean() {
        return root(PI * sigma / 2, 2);
    }

    double Median() {
        return root(sigma * 2 * ln(2), 2);
    }

    double Mode() {
        return sigma;
    }

    double Variance() {
        return ((4 - PI) / 2) * sigma;
    }

    double Skewness() {
        return (2 * root(PI, 2) * (PI - 3)) / ((4 - PI) * root(4 - PI, 2));
    }

    double EXkurtosis() {
        double part = -1 * (6 * power(PI, 2) - (24 * PI) + 16);
        return part / (power(4 - PI, 2));
    }

    double Entropy() {
        return 1 + ln(root(sigma / 2, 2)) + (EulerMascheroni / 2);
    }

    double MGF(double t) {
        double part1 = 1 + (root(sigma, 2) * t * euler(sigma * power(t, 2) / 2));
        double part2 = root(PI / 2, 2) * (errorFunction(root(sigma / 2, 2) * t) + 1);
        return part1 * part2;
    } 
}
