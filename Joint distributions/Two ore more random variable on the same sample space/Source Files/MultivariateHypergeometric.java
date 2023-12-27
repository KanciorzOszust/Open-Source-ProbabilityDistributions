public class MultivariateHypergeometric extends MathLibrary{
    int c;
    int[] K;
    int N;
    int n;

    MultivariateHypergeometric(int c, int[] K, int n) {
        if (c <= 0) throw new IllegalArgumentException("c > 0");
        if (K.length != c) throw new IllegalArgumentException("k.lenght = c");
        int N = 0;
        for (int i = 0; i < c; i++) {
            N += K[i];
        }
        if (n < 0 || n > N) throw new IllegalArgumentException("0 <= n <= sum of K");
        this.c = c;
        this.K = K;
        this.N = N;
        this.n = n;
    }

    double PMF(int[] k) {
        if (k.length != c) throw new IllegalArgumentException("k.length = c");
        double value = 1;
        for (int i = 0; i < c; i++) {
            value *= binomial(K[i], k[i]);
        }
        return value / binomial(N, n);
    }

    double Entropy(int i) {
        return n * (K[i] / N);
    }

    double Variance(int i) {
        double part = ((N - n) / (N - 1)) * (K[i] / N);
        return n * part * (1 - (K[i] / N));
    }
}
