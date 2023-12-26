public class TukeyLambda extends MathLibrary{
    double lambda;

    TukeyLambda(double lambda) {
        this.lambda = lambda;
    }

    double Quantile(double p) {
        if (lambda != 0) {
            return (1 / lambda) * (doublePower(p, lambda) - doublePower(1 - p, lambda));
        } else {
            return ln(p / (1 - p));
        }
    }

    double CDF(double x) {
        if (lambda == 0) {
            return 1 / (euler(-1 * x) + 1);
        } else {
            return Double.NaN;
        }
    }

    double Mean() {
        if (lambda > -1) return 0;
        else return Double.NaN;
    }

    int Median() {
        return 0;
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        if (lambda == 0) {
            return power(PI, 2) / 3;
        } else if (lambda > -0.5) {
            double part = (1 / (1 + (2 * lambda))) - (power(gamma(lambda + 1), 2) / gamma(2 * lambda + 2));
            return (2 / power(lambda, 2)) * part;
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (lambda > -1.0 / 3) return 0;
        else return Double.NaN;
    }

    private double g(int k) {
        return gamma(k * lambda + 1);
    }

    double Kurtosis() {
        if (lambda == 0) {
            return 1.2;
        } else if (lambda > -0.25) {
            double part1 = power(2 * lambda, 2) / (2 * (4 * lambda + 1));
            double part2 = power(g(2), 2) * (power(g(0), 2) - (4 * g(1) * g(3)) + g(4));
            double part3 = g(4) * power(power(g(1), 2) - g(2), 2);
            return part1 * (part2 / part3) - 3;
        } else {
            return Double.NaN;
        }
    }
}
