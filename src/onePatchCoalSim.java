import cern.jet.random.Binomial;
import cern.jet.random.Exponential;
import cern.jet.random.Uniform;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.MathUtils;
import sun.awt.datatransfer.DataTransferer;

import java.lang.reflect.Array;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/*

adapted from phydyn R package from Karcher et al - https://github.com/mdkarcher/phylodyn/tree/0750b0fd74ce44f42838f7dfd52df9d4f7fa5e6e/R
mainly from coalsieve.R and calculate.R


This is a test class that aims to simulate a coalescent tree under logistic growth population

 */

public class onePatchCoalSim {

    int K;
    double r;
    int N_init;

    public onePatchCoalSim() {

        K = 500;
        r = 0.05;
        N_init = 1;

    }

    public double Nt(double t) {

        double denominator = K + N_init * Math.exp(r * t) - 1.0; // K * P_0(e(-rt) - 1)
        double numerator = K * N_init * Math.exp(r * t); // K * P_0 * e(rt)

        return numerator / denominator;
    }

    /*
    sampleTimes = the times at which lineages are are sampled
    nSampled = the number of lineages at each sampled timepoint
    lowerBound = the minimum size of the population
     */

    public coalData coalSim_thin(List<Double> sampleTimes, List<Integer> nSampled, Double lowerBound) {

        List<Double> coalTimes = new ArrayList<>();
        List<Integer> lineages = new ArrayList<>();
        List<Double> intercoalTimes = new ArrayList<>();

        int curr = 0;
        int active_lineages = nSampled.get(curr);
        double time = sampleTimes.get(curr);

        while (time <= Collections.max(sampleTimes) || active_lineages > 1) {

            if (active_lineages == 1) {
                curr = curr + 1;
                active_lineages = active_lineages + nSampled.get(curr);
                time = sampleTimes.get(curr);

            }

            time = time + Exponential.staticNextDouble(0.5 * active_lineages * (active_lineages - 1) / lowerBound);

            if (curr < sampleTimes.size() - 1 && time >= sampleTimes.get(curr + 1)) {

                curr = curr + 1;
                active_lineages = active_lineages + nSampled.get(curr);
                time = sampleTimes.get(curr);

            } else if (Uniform.staticNextDouble() <= lowerBound / Nt(time)) {

                coalTimes.add(time);
                lineages.add(active_lineages);
                active_lineages = active_lineages - 1;

            }
        }

        coalData coalData = new coalData();
        coalData.setCoalTimes(coalTimes);
        coalData.setLineages(lineages);
        coalData.setIntercoalTimes(diff(coalTimes));
        coalData.setSampleTimes(sampleTimes);
        coalData.setnSampled(nSampled);

        return coalData;
    }

    public List<Double> diff(List<Double> coalTimes) {

        List<Double> intercoalTimes = new ArrayList<>();
        intercoalTimes.add(coalTimes.get(0));
        for (int i = 0; i < coalTimes.size() - 1; i++) {

            intercoalTimes.add(coalTimes.get(i + 1) - coalTimes.get(i));

        }

        return intercoalTimes;
    }

    public args2 gen_INLA_args(List<Double> sampleTimes, List<Integer> nSampled, List<Double> coalTimes) {

        args2 args2 = new args2();

        if (nSampled.stream().reduce(0, Integer::sum) != coalTimes.size() + 1) {

            throw new ArithmeticException("Number sampled not equal to number of coalescent events + 1.");

        }

        if (coalTimes.stream()
                .distinct()
                .filter(sampleTimes::contains)
                .collect(Collectors.toSet()).size() > 0) {


            System.out.println("Coincident sampling event and coalescent event: results may be unpredictable.");
        }

        int l = sampleTimes.size();
        int m = coalTimes.size();
        List<Double> sorting = Stream.of(sampleTimes, coalTimes).flatMap(x -> x.stream()).collect(Collectors.toList());


        List<Integer> indices = sortIndicesByAnotherList(sorting);

        //lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix]

        List<Integer> lineage_change = sortListByIndices(indices, Stream.of(nSampled, rep(m, -1)).
                flatMap(x -> x.stream()).
                collect(Collectors.toList()));



        List<Integer> lineages = IntStream.range(0, lineage_change.size() - 1)
                .map(i -> IntStream.rangeClosed(0, i).map(lineage_change::get).sum())
                .boxed()
                .collect(Collectors.toList());

        List<Double> coalFactor = lineages.stream().map(i -> (double)i * (i - 1) / 2.0).collect(Collectors.toList());

        //System.out.println(Arrays.asList(Stream.of(rep(l, 0), rep(m, 1)).toArray()));
        List<Integer> event = sortListByIndices(indices, Stream.of(rep(l, 0), rep(m, 1)).
                flatMap(x -> x.stream()).
                collect(Collectors.toList()));



        args2.setCoalFactor(coalFactor);
        args2.setSorting(sorting);
        args2.setEvent(event);
        args2.setLineages(lineages);

        return args2;

    }

    public List<Integer> rep(int n, int value) {


        List<Integer> a = new ArrayList<>(Arrays.asList(new Integer[n]));
        Collections.fill(a,value);

        return a;
    }

    public List<Integer> sortIndicesByAnotherList(List<Double> list) {

        List<Integer> indices = IntStream.range(0, list.size()).
                boxed().
                collect(Collectors.toList());

        Collections.sort(indices, Comparator.comparingDouble(list::get));
        Collections.sort(list);

        return indices;
    }

    public List<Integer> sortListByIndices(List<Integer> indices, List<Integer> list) {

        List<Integer> sortedList = new ArrayList<>();
        indices.forEach(s -> sortedList.add(list.get(s)));

        return sortedList;

    }

    public void sampleTree(coalData gene) {

        int n = (gene.nSampled).stream().mapToInt(i -> i.intValue()).sum();
        int Nnode = n - 1;

        List<String> labels = new ArrayList<>();
        IntStream.range(0,n).forEach(i -> labels.add("t"+(i+1)));

        int tb = gene.nSampled.get(0); // Total branches (initial)
        double s = 0; //time for branch lengths;
        List<String> temp_labels = new ArrayList<>();
        temp_labels.addAll(labels.subList(0, tb));

        List<Double> temp_times = new ArrayList(Collections.nCopies(gene.nSampled.get(0), gene.sampleTimes.get(0)));

        int initial_row = 1;

        args2 args2 = gen_INLA_args(gene.sampleTimes, gene.nSampled, gene.coalTimes);

        System.out.println(gene.coalTimes);

        for (int j=1; j < args2.event.size(); j++) {

            //System.out.println("e "+j+", "+args2.event.get(j));
            if (args2.event.get(j) == 1) {

                s = args2.sorting.get(j);

                List<Integer> ra = sample(tb, 2);

                String new_label = "(" + temp_labels.get(ra.get(0)) + ":" + (s-temp_times.get(ra.get(0))) + ","
                        + temp_labels.get(ra.get(1)) + ":" + (s-temp_times.get(ra.get(1))) + ")";


                temp_labels.set(ra.get(0), new_label);
                temp_labels.remove((int)ra.get(1));

                temp_times.set(ra.get(0), s);
                temp_times.remove((int) ra.get(1));
                tb = tb - 1;
                System.out.println(temp_labels.get(ra.get(0)));

            }
            else {
                //I will be adding samples at
                s = args2.sorting.get(j);

                if (gene.nSampled.get(initial_row) == 1) {

                    temp_labels.add(labels.get(cumsum(gene.nSampled).get(initial_row)-1));
                    initial_row = initial_row + 1;
                    tb = tb + 1;
                    temp_times.add(s);
                }
                else {
                    int end = cumsum(gene.nSampled).get(initial_row);
                    int ini = cumsum(gene.nSampled).get(initial_row - 1);
                    for (int k = ini; k < end; k++) {
                        temp_labels.add(labels.get(k));
                        tb = tb + 1;
                        temp_times.add(s);
                    }
                    initial_row = initial_row + 1;
                }
            }
        }
    }

    public static void main (String [] args) {

        onePatchCoalSim a = new onePatchCoalSim();


        List<Double> sampleTimes = new ArrayList<>(Arrays.asList(0.0,10.0,20.0));
        List<Integer> nSampled = new ArrayList<>(Arrays.asList(10,10,10));
        double lowerBound = 1.0;

        coalData gene = a.coalSim_thin(sampleTimes, nSampled, lowerBound);

        a.sampleTree(gene);


    }

    //

    private List<Integer> cumsum(List<Integer> list) {

        List<Integer> cumsum = IntStream.range(0, list.size())
                .map(i -> IntStream.rangeClosed(0, i).map(list::get).sum())
                .boxed()
                .collect(Collectors.toList());

        return cumsum;

    }

    private List<Integer> sample(int tb, int n) {

        List<Integer> indices = IntStream.range(0, tb).
                boxed().
                collect(Collectors.toList());

        Collections.shuffle(indices);

        List<Integer> temp = new ArrayList<>();

        temp.addAll(indices.subList(0, n));
        Collections.sort(temp);
        return temp;

    }

    class args2 {

        List<Double> coalFactor = new ArrayList<>();
        List<Double> sorting = new ArrayList<>();
        List<Integer> event = new ArrayList<>();
        List<Integer> lineages = new ArrayList<>();

        void setCoalFactor(List<Double> coalFactor) {

            this.coalFactor = coalFactor;
        }
        void setSorting(List<Double> sorting) {

            //Collections.reverse(sorting);
            System.out.println(sorting);
            this.sorting = sorting;
        }
        void setEvent(List<Integer> event) {

            this.event = event;
        }
        void setLineages(List<Integer> lineages) {

            this.lineages = lineages;
        }

        List<Double> getCoalFactor() {

            return this.coalFactor;
        }
        List<Double> getSorting() {

            return this.sorting;
        }
        List<Integer> getEvent() {

            return this.event;
        }
        List<Integer> getLineages() {

            return this.lineages;
        }



    }

    class coalData {

        List<Double> coalTimes = new ArrayList<>();
        List<Integer> lineages = new ArrayList<>();
        List<Double> intercoalTimes = new ArrayList<>();
        List<Double> sampleTimes = new ArrayList<>();
        List<Integer> nSampled = new ArrayList<>();


        public List<Double> getCoalTimes() {
            return coalTimes;
        }

        public void setCoalTimes(List<Double> coalTimes) {
            this.coalTimes = coalTimes;
        }

        public List<Integer> getLineages() {
            return lineages;
        }

        public void setLineages(List<Integer> lineages) {
            this.lineages = lineages;
        }

        public List<Double> getIntercoalTimes() {
            return intercoalTimes;
        }

        public void setIntercoalTimes(List<Double> intercoalTimes) {
            this.intercoalTimes = intercoalTimes;
        }

        public List<Double> getSampleTimes() {
            return sampleTimes;
        }

        public void setSampleTimes(List<Double> sampleTimes) {
            this.sampleTimes = sampleTimes;
        }

        public List<Integer> getnSampled() {
            return nSampled;
        }

        public void setnSampled(List<Integer> nSampled) {
            this.nSampled = nSampled;
        }
    }
}
