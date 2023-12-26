public class InverseGamma extends MathLibrary{
    double alpha;
    double beta;

    InverseGamma(double alpha, double beta) {
        if (alpha <= 0 || beta <= 0) throw new IllegalArgumentException("alpha, beta > 0");
        this.alpha = alpha;
        this.beta = beta;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part = doublePower(beta, alpha) / gamma(alpha) * doublePower(x, -1 * alpha - 1);
        return part * euler(-1 * beta / x);
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return incompleteUpperGamma(alpha, beta / x) / gamma(alpha);
    }

    double Mean() {
        if (alpha > 1) return beta / (alpha - 1);
        else return Double.NaN;
    }

    double Mode() {
        return beta / (alpha + 1);
    }

    double Variance() {
        if (alpha > 2) {
            return power(beta, 2) / (power(alpha - 1, 2) * (alpha - 2));
        } else {
            return Double.NaN;
        }
    }

    double Skewness() {
        if (alpha > 3) {
            return (4 * root(alpha - 2, 2)) / (alpha - 3);
        } else {
            return Double.NaN;
        }
    }

    double EXkurotsis() {
        if (alpha > 4) {
            return (6 * (5 * alpha - 11)) / ((alpha - 3) * (alpha - 4));
        } else {
            return Double.NaN;
        }
    }

    double Entropy() {
        return alpha + (beta * gamma(alpha)) - ((1 + alpha) * digamma(alpha));
    }
} 
