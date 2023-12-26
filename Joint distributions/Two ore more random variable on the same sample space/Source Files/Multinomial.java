public class Multinomial extends MathLibrary{
    int n;
    int k;
    double[] p;

    Multinomial(int n, double[] p) {
        if (n <= 0) throw new IllegalArgumentException("n > 0");
        if (p.length == 0) throw new IllegalArgumentException("p.length > 0");
        double value = 0;
        for (int i = 0; i < p.length; i++) {
            value += p[i];
        } 
        if (value != 1) throw new IllegalArgumentException("sum of p[i] = 1");
        this.n = n;
        this.k = p.length;
        this.p = p;
    }

    double PMF(double[] x) {
        if (x.length != k) throw new IllegalArgumentException("x.length = p.length");
        double value = 0;
        for (int i = 0; i < k; i++) {
            if (x[i] < 0) throw new IllegalArgumentException("x[i] >= 0");
            value += x[i];
        }
        if (value != n) throw new IllegalArgumentException("sum of x[i] == n");

        double value1 = 1;
        double value2 = 1;
        for (int i = 0; i < k; i++) {
            value *= gamma(x[i] + 1);
        }
        for (int i = 0; i < k; i++) {
            value2 *= doublePower(p[i], x[i]);
        }

        return factorial(n) / value1 * value2;
    }

    double Mean(int i) {
        return n * p[i];
    }

    double Variance(int i) {
        return Mean(i) * (1 - p[i]);
    }

    double Entropy(int i) {
        double value1 = 0;
        double value2 = 0;
        for (int j = 0; j < k; j++) {
            value1 += p[k] * ln(p[k]);
        }
        for (int j = 0; j < k; j++) {
            for (int j2 = 0; j2 < n; j2++) {
                value2 += binomial(n, j2) * power(p[i], j2) * power(1 - p[i], n - j2) * ln(factorial(j2));
            }
        }
        return (-1 * ln(factorial(n))) - (n * value1) + value2;
    }

    double MGF(double[] t) {
        double value = 0;
        for (int i = 0; i < k; i++) {
            value += p[i] * euler(t[i]);
        }
        return power(value, n);
    }
}
