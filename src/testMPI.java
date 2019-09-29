import mpi.*;


import java.nio.ByteBuffer;


public class testMPI {

    public testMPI() {

    }

    public static void main(String [] args) throws MPIException {

        int rank, size, tag;
        int start,stop,i;
        int left, right, sum;
        ByteBuffer bb;
        int[] passon = new int[1];
        int[] addon = new int[1];

        Status status;
        Request request;


        tag = 1;

        MPI.Init(args);
        rank = MPI.COMM_WORLD.getRank();
        size = MPI.COMM_WORLD.getSize();

        //Send current number to the right and receive from the left
        left = rank -1;
        right = rank + 1;
        if ((rank+1) == size){
            right = 0;
        }
        if(rank == 0){
            left = size - 1;
        }

        sum = 0;

        // Initialise local values to:
        passon[0] = ((rank+1)*(rank+1));

        //  Alternatively use
        //passon[0] = rank;

        // Use non-blocking point-to-point communication

        for(i=0;i<size;i++){

            //if(rank == i) {
                MPI.COMM_WORLD.sendRecv(passon, passon.length, MPI.INT, right, 1, addon, addon.length, MPI.INT, left, 1);
                sum = sum + addon[0];
                System.out.println(sum+", "+addon[0]);
                passon = addon;
            //}

        }


        System.out.printf("The sum is: %d on processor %d \n",sum, rank);

        MPI.Finalize();

    }

//    public static void main(String [] args) throws MPIException {
//        double pi = 0, exactpi = 0.0;
//
//        int i;
//
//        /* MPI variables */
//
//
//        int rank, size, source, tag;
//
//        /* Other variables */
//
//        int istart, istop;
//        double[] partialpi = new double[1], recvpi = new double[1];
//
//
//        /* Initialise MPI and compute number of processes and local rank */
//
//
//        MPI.Init(args);
//        size = MPI.COMM_WORLD.getSize();
//        rank = MPI.COMM_WORLD.getRank();
//
//        int N = 10;
//        if (rank == 0) {
//            System.out.printf("\nComputing approximation to pi using N = %d", N);
//        }
//
//        /* Now make sure output only comes from one process */
//
//        if (rank == 0) {
//            System.out.printf("\nRunning on %d process(es)", size);
//        }
//
//        partialpi[0] = 0.0;
//
//        /*
//         * Compute an approximation to pi using a simple series expansion for pi/4
//         * Ensure each process computes a separate section of the summation
//         * NOTE: here I assume that N is exactly divisible by the number of processes
//         */
//
//        ArrayList<Double> b = new ArrayList<>();
//
//        for(int t = 0; t < 2; t++) {
//
//            if (rank == 0) {
//                partialpi[0] += 2;
//
//            }
//
//            b.add(1.0);
//            b.add(1.0);
//
//
//            for (i = 0; i < b.size(); i++) {
//
//
//                if (rank == 0) {
//
//                    for (source = 1; source < 2; source++) {
//                        /* receive partialpi from rank=source and place value in recvpi */
//                        /* all messages are tagged as zero */
//
//                        tag = 0;
//
//                        MPI.COMM_WORLD.recv(recvpi, 1, MPI.DOUBLE, source, tag);
//
//                        /* add to running total */
//
////                System.out.println("\n>>>"+rank+", "+recvpi[0]);
//                        pi = recvpi[0];
//
//                        //System.out.println(i + ", " + pi+", "+rank);
//                    }
//                }
//                else if(rank == 1) {
//                    partialpi[0] += 1;
//                    tag = 0;
//                    MPI.COMM_WORLD.sSend(partialpi, 1, MPI.DOUBLE, 0, tag);
//                    //System.out.println(i + " " + partialpi[0] + " " + rank);
//                }
//
//            }
//            b.clear();
//
//            istart = (N/size) * rank + 1;
//            istop = istart + (N/size) - 1;
//
//            for (i = istart; i <= istop; i++) {
//
//                partialpi[0] += 1;
//
//            }
//            if (rank == 0) {
//
//                pi = partialpi[0];
//                for (source = 1; source < size; source++) {
//                    /* receive partialpi from rank=source and place value in recvpi */
//                    /* all messages are tagged as zero */
//
//                    tag = 0;
//
//                    //System.out.println(istart+", "+istop+", "+recvpi[0]+", "+rank);
//                    MPI.COMM_WORLD.recv(recvpi, 1, MPI.DOUBLE, source, tag);
//                    System.out.println(N+", "+recvpi[0]+", "+rank);
//
//                    /* add to running total */
//
////                System.out.println("\n>>>"+rank+", "+recvpi[0]);
//                    pi += recvpi[0];
//
//                    //System.out.println(i + ", " + pi+", "+rank);
//                }
//            }
//            else {
//
//                tag = 0;
//                MPI.COMM_WORLD.sSend(partialpi, 1, MPI.DOUBLE, 0, tag);
//                //System.out.println(i + " " + partialpi[0] + " " + rank);
//            }
//
//
//        }
//        System.out.printf("\nOn rank %d partialpi = %f, and pi = %f", rank, partialpi[0], pi);
//
//        /*
//         * Compute global value of pi by sending partial values to rank 0
//         * NOTE: this would be more efficiently done using MPI_REDUCE
//         */
//
////        if (rank == 0) {
////            /* Initialise pi to locally computed parial sum */
////
////            System.out.println("\n>"+rank+", "+partialpi[0]);
////            pi = partialpi[0];
////
////            /* Add in contribution from other processes */
////
////            for (source = 1; source < size; source++) {
////                /* receive partialpi from rank=source and place value in recvpi */
////                /* all messages are tagged as zero */
////
////                tag = 0;
////
////
////                MPI.COMM_WORLD.recv(recvpi, 1, MPI.DOUBLE, source, tag);
////
////                /* add to running total */
////
//////                System.out.println("\n>>>"+rank+", "+recvpi[0]);
////                pi = pi + recvpi[0];
////            }
////        }
////        else {
////            /* all other processes send their partial value to rank 0 */
////
////            tag = 0;
////
////            MPI.COMM_WORLD.sSend(partialpi, 1, MPI.DOUBLE, 0, tag);
////        }
//
////        pi = (pi * 4.0) / ((double) N);
////
////        exactpi = 4.0 * Math.atan(1.0);
////
////        if (rank == 0) {
////            System.out.println("\npi = " + pi + " % error = " + Math.abs((100.0 * (pi - exactpi))/exactpi));
////        }
//
//        MPI.Finalize();
//
//
//    }
}


