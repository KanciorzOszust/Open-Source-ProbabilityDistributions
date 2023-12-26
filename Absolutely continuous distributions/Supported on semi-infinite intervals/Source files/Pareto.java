public class Pareto extends MathLibrary{
    double xm;
    double alpha;

    Pareto(double xm, double alpha) {
        if (xm <= 0 || alpha <= 0) throw new IllegalArgumentException("xm, alpha > 0");
        this.xm = xm;
        this.alpha = alpha;
    }

    double PDF(double x) {
        if (x < xm) throw new IllegalArgumentException("x >= xm");
        return (alpha * doublePower(xm, alpha)) / doublePower(x, alpha + 1);
    }

    double CDF(double x) {
        if (x < xm) throw new IllegalArgumentException("x >= xm");
        return 1 - doublePower(xm / x, alpha);
    }

    double Quantile(double p) {
        return xm * doublePower(1 - p, -1 / alpha);
    }

    double Mean() {
        if (alpha > 1) {
            return (alpha * xm) / (alpha - 1);
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Median() {
        return xm * doublePower(2, 1 / alpha);
    }

    double Mode() {
        return xm;
    }

    double Variance() {
        if (alpha > 2) {
            return (power(xm, 2) * alpha) / (power(alpha - 1, 2) * (alpha - 2));
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Skewness() {
        if (alpha > 3) {
            return ((2 * (1 + alpha)) / (alpha - 3)) * root((alpha - 2) / alpha, 2);
        } else {
            return Double.NaN;
        }
    }

    double EXkurtosis() {
        if (alpha > 4) {
            double part1 = 6 * (power(alpha, 3) + power(alpha, 2) - (6 * alpha) - 2);
            double part2 = alpha * (alpha - 3) * (alpha - 4);
            return part1 / part2;
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        return log((xm / alpha) * euler(1 + (1 / alpha)), 10);
    }

    double[][] FisherInformation() {
        double[][] values = {
            {power(alpha, 2) / power(xm, 2), 0},
            {0, 1 / power(alpha, 2)}
        };
        return values;
    }

    double ExpectedShortfall(double p) {
        return (xm * alpha) / (doublePower(1 - p, 1 / alpha) * (alpha - 1));
    }
}
