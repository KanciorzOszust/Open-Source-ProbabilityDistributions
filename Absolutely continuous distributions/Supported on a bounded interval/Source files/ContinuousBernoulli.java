public class ContinuousBernoulli extends MathLibrary{
    double lambda;

    ContinuousBernoulli(double lambda) {
        if (lambda <= 0 || lambda >= 1) throw new IllegalArgumentException(" 0 < lambda < 1");
        this.lambda = lambda;
    }

    private double C(double x) {
        if (x== 0.5) {
            return 2;
        } else {
            return (2 * hyperbolicArcTangent(1 - (2 * lambda))) / (1 - (2 * lambda));
        }
    }

    double PDF(double x) {
        if (x < 0 || x > 1) throw new IllegalArgumentException("0 <= x <= 1");
        return C(lambda) * doublePower(lambda, x) * doublePower(1 - lambda, 1 - x);
    }

    double CDF(double x) {
        if (lambda == 0.5) {
            return x;
        } else {
            double part = doublePower(lambda, x) * doublePower(1 - lambda, 1 - x) + lambda - 1;
            return part / (2 * lambda - 1);
        }
    }

    double Mean() {
        if (lambda == 0.5) {
            return 0.5;
        } else {
            return (lambda / (2 * lambda - 1)) + (1 / (2 * hyperbolicArcTangent(1 - (2 * lambda))));
        }
    }

    double Variance() {
        if (lambda == 0.5) {
            return 1.0 / 12;
        } else {
            double part1 = (-1 * (1 - lambda) * lambda) / power(1 - (2 * lambda), 2);
            double part2 = 1 / power(2 * hyperbolicArcTangent(1 - (2 * lambda)), 2);
            return part1 + part2;
        }
    }
}
