public class StudentsT extends MathLibrary{
    double nu;

    StudentsT(double nu) {
        if (nu <= 0) throw new IllegalArgumentException("nu > 0");
        this.nu = nu;
    }

    double PDF(double x) {
        double part1 = gamma((nu + 1) / 2) / (root(nu * PI, 2) * gamma(nu / 2));
        double part2 = doublePower(1 + (power(x, 2) / nu), -1 * (nu + 1) / 2);
        return part1 * part2;
    }

    double CDF(double x) {;
        double part1 = x * gamma((nu + 1) / 2);
        double part2 = hypergeometric(0.5, (nu + 1) / 2, 1.5, -1 * power(x, 2) / nu) / (root(PI * nu, 2) * gamma(nu / 2));
        return part1 * part2;
    }

    double Mean() {
        if (nu > 1) return 0;
        else return Double.NaN;
    }

    int Median() {
        return 0;
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        if (nu > 2) return nu / (nu - 2);
        else if (nu > 1) return Double.POSITIVE_INFINITY;
        else return Double.NaN;
    }

    double Skewness() {
        if (nu > 3) return 0;
        else return Double.NaN;
    }

    double EXkurtosis() {
        if (nu > 4) return 6 / (nu - 4);
        else if (nu > 2) return Double.POSITIVE_INFINITY;
        else return Double.NaN;
    }

    double Entropy() {
        double part1 = ((nu +1) / 2) * (digamma((1 + nu) / 2) - digamma(nu / 2));
        double part2 = ln(root(nu, 2) * beta(nu / 2, 0.5));
        return part1 + part2;
    }
}
