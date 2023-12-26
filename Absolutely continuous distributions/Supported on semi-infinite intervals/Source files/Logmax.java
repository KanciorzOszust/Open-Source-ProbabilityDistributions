public class Logmax extends MathLibrary{
    double alpha;
    double gamma;

    Logmax(double alpha, double gamma) {
        if (alpha <= 0 || gamma <= 0) throw new IllegalArgumentException("alpha, gamma > 0");
        this.alpha = alpha;
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (alpha / gamma) * doublePower(1 + (x / gamma), -1 * (alpha + 1)); 
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - doublePower(1 + (x / gamma), -1 * alpha);
    }

    double Quantile(double p) {
        return gamma * (doublePower(1 - p, -1 / alpha) - 1);
    }

    double Mean() {
        if (alpha > 1) return gamma / (alpha - 1);
        else return Double.NaN;
    }

    double Median() {
        return gamma * (doublePower(2, 1 / alpha) - 1);
    }

    int Mode() {
        return 0;
    }

    double Variance() {
        if (alpha > 2) {
            return (power(gamma, 2) * alpha) / (power(alpha - 1, 2) * (alpha - 2));
        } else if (alpha <= 2 && alpha > 1) {
            return Double.POSITIVE_INFINITY;
        } else {
            return Double.NaN;
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
}
