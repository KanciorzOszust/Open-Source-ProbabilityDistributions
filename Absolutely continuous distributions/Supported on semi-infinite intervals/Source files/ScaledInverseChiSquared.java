public class ScaledInverseChiSquared extends MathLibrary{
    double nu;
    double r;

    ScaledInverseChiSquared(double nu, double r) {
        if (nu <= 0 || r <= 0) throw new IllegalArgumentException("nu, r > 0");
        this.nu = nu;
        this.r = r;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = doublePower(r * nu / 2, nu / 2) / gamma(nu / 2);
        double part2 = euler(-1 * (nu * r) / (2 * x)) / doublePower(x, 1 + (nu / 2));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return incompleteUpperGamma(nu / 2, (r * nu) / (2 * x)) / gamma(nu / 2);
    }

    double Mean() {
        if (nu > 2) return (nu * r) / (nu - 2);
        else return Double.NaN;
    }

    double Mode() {
        return (nu * r) / (nu + 2);
    }

    double Variance() {
        if (nu > 4) {
            return (2 * power(nu, 2) * power(r, 2)) / (power(nu - 2, 2) * (nu - 8));
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
        double part1 = (nu / 2) + ln((r * nu) / 2 * gamma(nu / 2));
        double part2 = (1 + (nu / 2)) * digamma(nu / 2);
        return part1 - part2;
    }

    double MGF(double t) {
        double part1 = 2 / gamma(nu / 2);
        double part2 = doublePower((-1 * r * nu * t) / 2, nu / 4);
        double part3 = modifiedBessel2(nu / 2, root(-2 * r * nu * t, 2));
        return part1 * part2 * part3;
    }
}