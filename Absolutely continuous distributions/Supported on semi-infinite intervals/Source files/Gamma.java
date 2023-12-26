public class Gamma extends MathLibrary{
    double k;
    double theta;

    Gamma(double k, double theta) {
        if (k <= 0 || theta <= 0) throw new IllegalArgumentException("k, theta > 0");
        this.k = k;
        this.theta = theta;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = 1 / (gamma(k) * doublePower(theta, k));
        double part2 = doublePower(x, k - 1) * euler(-1 * x / theta);
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return incompleteLowerGamma(k, x / theta) / gamma(k);
    }

    double Mean() {
        return k * theta;
    }

    double Mode() {
        if (k >= 1) return (k - 1) * theta;
        else return 0;
    }

    double Variance() {
        return k * power(theta, 2);
    }

    double Skewness() {
        return 2 / root(k, 2);
    }

    double EXkurtosis() {
        return 6 / k;
    }

    double Entropy() {
        return k + ln(theta) + ln(gamma(k)) + ((1 - k) * digamma(k));
    }

    double MGF(double t) {
        if (t < 1 / theta) {
            return doublePower(1 - (theta * k), -1 * k);
        } else {
            return Double.NaN;
        }
    }

    double[][] FIsherInformation() {
        double[][] values = {
            {trigamma(k), 1 / theta},
            {1 / theta, k / power(theta, 2)}
        };
        return values;
    }
}
