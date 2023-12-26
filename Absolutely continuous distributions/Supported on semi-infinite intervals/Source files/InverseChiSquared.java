public class InverseChiSquared extends MathLibrary{
    double nu;

    InverseChiSquared(double nu) {
        if (nu <= 0) throw new IllegalArgumentException("nu > 0");
        this.nu = nu;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = doublePower(2, -1 * nu / 2) / gamma(nu / 2);
        double part2 = doublePower(x, (-1 * nu / 2) - 1) * euler(-1 / (2 * x));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return incompleteUpperGamma(nu / 2, 1 / (2 * x)) / gamma(nu / 2);
    }

    double Mean() {
        if (nu > 2) return 1 / (nu - 2);
        else return Double.NaN;
    }

    double Median() {
        return 1 / (nu * power(1 - (2 / (9 * nu)), 3));
    }

    double Mode() {
        return 1 / (nu + 2);
    }

    double Variance() {
        if (nu > 4) {
            return 2 / (power(nu - 2, 2) * (nu - 4));
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (nu > 6) {
            return (4 / (nu - 6)) * root(2 * (nu - 4), 2);
        } else {
            return Double.NaN;
        }
    }

    double EXkurtosis() {
        if (nu > 8) {
            return (15 * (5 * nu - 22)) / ((nu - 6) * (nu - 8));
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        double part1 = (nu / 2) + ln((nu / 2) * gamma(nu / 2));
        double part2 = (1 + (nu / 2)) * digamma(nu / 2);
        return part1 - part2;
    }
}
