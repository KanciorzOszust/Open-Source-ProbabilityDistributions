public class NoncentralChiSquared extends MathLibrary{
    double k;
    double gamma;

    NoncentralChiSquared(double k, double gamma) {
        if (k <= 0 || gamma <= 0) throw new IllegalArgumentException("k, gamma > 0");
        this.k = k;
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = 0.5 * euler(-1 * (x * gamma) / 2);
        double part2 = doublePower(x / gamma, (k / 4) - 0.5);
        double part3 = modifiedBessel((k / 2) - 1, root(gamma * x, 2));
        return part1 * part2 * part3;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - marcumQ(k / 2, root(gamma, 2), root(x, 2));
    }

    double Mean() {
        return k + gamma;
    }

    double Variance() {
        return 2 * (k + (2 * gamma));
    }

    double Skewness() {
        double part1 = 2 * root(2, 2) * (k + (3 * gamma));
        double part2 = (k + (2 * gamma)) * root(k + (2 * gamma), 2);
        return part1 / part2;
    }

    double EXkurtosis() {
        return (12 * (k + (4 * gamma))) / power(k + (2 * gamma), 2);
    }

    double MGF(double t) {
        if (t < 0.5) {
            return euler((gamma * t) / (1 - (2 * t))) / doublePower(1 - (2 * t), k /2);
        } else {
            return Double.NaN;
        }
    }
}
