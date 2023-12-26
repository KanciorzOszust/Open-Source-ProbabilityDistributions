public class GeneralizedPareto extends MathLibrary{
    double mu;
    double sigma;
    double eta;

    GeneralizedPareto(double mu, double sigma, double eta) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
        this.eta = eta;
    }

    private double z(double x) {
        return (x - mu) / sigma;
    }

    double PDF(double x) {
        if (eta >= 0) {
            if (x < mu) throw new IllegalArgumentException("x >= mu");
        } else {
            if (x < mu || x > (mu - (sigma / eta))) throw new IllegalArgumentException("x >= mu");
        }
        return (1 / sigma) * doublePower(1 + (eta * z(x)), -1 * (1 / (eta + 1)));
    }

    double CDF(double x) {
        if (eta >= 0) {
            if (x < mu) throw new IllegalArgumentException("x >= mu");
        } else {
            if (x < mu || x > (mu - (sigma / eta))) throw new IllegalArgumentException("x >= mu");
        }
        return 1 - doublePower(1 + (eta * z(x)), -1 / eta);
    }

    double Mean() {
        if (eta < 1) {
            return mu + (sigma / (1 - eta));
        } else {
            return Double.NaN;
        }
    }

    double Median() {
        return mu + ((sigma * (doublePower(2, eta) - 1)) / eta);
    }

    double Mode() {
        return eta;
    }

    double Variance() {
        if (eta < 0.5) {
            return power(sigma, 2) / (power(1 - eta, 2) * (1 - (2 * eta)));
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (eta < 1.0 / 3) {
            return (2 * (1 + eta) * root(1 - (2 * eta), 2)) / (1 - (3 * eta));
        } else {
            return Double.NaN;
        }
    }

    double EXkurtosis() {
        if (eta < 0.25) {
            double part1 = 3 * (1 - (2 * eta)) * (2 * power(eta, 2) + eta + 3);
            double part2 = (1 - (3 * eta)) * (1 - (4 * eta));
            return part1 / part2;
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        return log(sigma, 10) + eta + 1;
    }

    double ExpectedShortfall(double p) {
        if (eta != 0) {
            double part1 = doublePower(1 - p, -1 * eta) / (1 - eta);
            double part2 = (doublePower(1 - p, -1 * eta) - 1) / eta;
            return mu + (sigma * (part1 + part2));
        } else {
            return mu + (sigma * (1 - ln(1 - p)));
        }
    }
}
