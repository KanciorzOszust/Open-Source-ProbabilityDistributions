public class Weibull extends MathLibrary{
    double lambda;
    double k;

    Weibull(double lambda, double k) {
        if (lambda <= 0 || k <= 0) throw new IllegalArgumentException("lambda, k > 0");
        this.lambda = lambda;
        this.k = k;
    }

    double PDF(double x) {
        if (x >= 0) {
            double part = (k / lambda) * doublePower(x / lambda, k - 1);
            return part * euler(-1 * doublePower(x / lambda, k));
        } else {
            return 0;
        }
    }

    double CDF(double x) {
        if (x >= 0) {
            return 1 - euler(-1 * doublePower(x / lambda, k));
        } else {
            return 0;
        }
    }

    double Mean() {
        return lambda * gamma(1 + (1/ k));
    }

    double Median() {
        return lambda * doublePower(ln(2), 1 / k);
    }

    double Mode() {
        if (k > 1) {
            return lambda * doublePower((k - 1) / k, 1 / k);
        } else {
            return 0;
        }
    }

    double Variance() {
        double part = gamma(1 + (2 / k)) - power(gamma(1 + (1 / k)), 2);
        return power(lambda, 2) * part;
    }

    double Skewness() {
        double part = power(lambda, 3) * gamma(1 + (3 / k)) - (3 * Mean() * Variance()) - power(Mean(), 3);
        return part / (Variance() * root(Variance(), 2));
    }

    double EXkurtosis() {
        double part1 = power(lambda, 4) * gamma(1 + (4 / k)) - (4 * Skewness() * Variance() * root(Variance(), 2) * Mean());
        double part2 = 6 * power(Mean(), 2) * Variance() - power(Mean(), 4);
        return (part1 - part2) / power(Variance(), 2) - 3;
    }

    double Entropy() {
        return EulerMascheroni * (1 - (1 / k)) + ln(lambda / k) + 1;
    }

    double MGF(double t) {
        if (k >= 1) {
            double value = 0;
            for (int n = 0; n < 7; n++) {
                value += power(t, n) * power(lambda, n) / factorial(n) * gamma((1 + n) / k);
            }
            return value;
        } else {
            return Double.NaN;
        }
    }
}
