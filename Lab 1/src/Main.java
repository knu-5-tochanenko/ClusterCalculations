import mpi.MPI;

public class Main {
    public static void main(String[] args) {
        MPI.Init(args);
        int me = MPI.COMM_WORLD.Rank();
        int size = MPI.COMM_WORLD.Size();
        System.out.println("Hello from " + me + ". Size: " + size);
        MPI.Finalize();
    }
}
