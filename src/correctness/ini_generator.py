class Args:
    def __init__(self, total_qbit = 12,mpi_qbit = 0 ,file_qbit = 6, chunk_qbit = 6,
            max_qbit = 38, max_path = 260, max_depth= 1000, is_density = 0,
            runner_type = "MEM", is_subcircuit = 1, is_directIO = 0, dump_file = "res.txt",MPI_testing = 0,MPI_buffe_size = 16,
            state_paths = "./state/path0,./state/path1,./state/path2,./state/path3,./state/path4,./state/path5,./state/path6,./state/path7,./state/path8,./state/path9,./state/path10,./state/path11,./state/path12,./state/path13,./state/path14,./state/path15,./state/path16,./state/path17,./state/path18,./state/path19,./state/path20,./state/path21,./state/path22,./state/path23,./state/path24,./state/path25,./state/path26,./state/path27,./state/path28,./state/path29,./state/path30,./state/path31,./state/path32,./state/path33,./state/path34,./state/path35,./state/path36,./state/path37,./state/path38,./state/path39,./state/path40,./state/path41,./state/path42,./state/path43,./state/path44,./state/path45,./state/path46,./state/path47,./state/path48,./state/path49,./state/path50,./state/path51,./state/path52,./state/path53,./state/path54,./state/path55,./state/path56,./state/path57,./state/path58,./state/path59,./state/path60,./state/path61,./state/path62,./state/path63"
        ):
        self.total_qbit = total_qbit
        self.mpi_qbit = mpi_qbit
        self.file_qbit = file_qbit
        self.chunk_qbit = chunk_qbit
        self.max_qbit = max_qbit
        self.max_path = max_path
        self.max_depth= max_depth
        self.is_density = is_density
        self.runner_type = runner_type
        self.is_subcircuit = is_subcircuit
        self.is_directIO = is_directIO
        self.dump_file = dump_file
        self.MPI_testing = MPI_testing
        self.MPI_buffe_size = MPI_buffe_size
        self.state_paths = state_paths


class Ini:
    def __init__(self, args, path) -> None:
        self.setting = args
        self.path = path

    def out(self):
        f=open(self.path,'w')
        print('[system]',file=f)
        args = self.setting.__dict__
        for key in args.keys():
            print(f'{key}={args[key]}',file=f)
        f.close()

# ************example************
# ini = Ini(Args(), "res.ini")
# ini.out()