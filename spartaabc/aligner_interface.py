import subprocess
    




class Aligner:
    def __init__(self, aligner: str) -> None:
        self._aligner_name = aligner

    def get_name(self):
        return self._aligner_name
    
    def set_input_file(self, file_path, tree_file=None):
        self._input_file = str(file_path)
        self._output_name = str(file_path).split(".")[0]
        if self._aligner_name == "MAFFT":
            self._aligner_cmd = ["mafftpy", "--strategy", "ginsi", self._input_file]
        if self._aligner_name == "MAFFTFAST":
            self._aligner_cmd = ["mafftpy", "--strategy", "fftns1" , self._input_file]


    def get_realigned_msa(self) -> str:
        # print(self._aligner_cmd)
        result = subprocess.run(self._aligner_cmd, capture_output=True, text=True)
        realigned_msa, stderr = result.stdout, result.stderr
        # realigned_msa = realigned_msa.split("\n")
        realigned_msa = [line for line in realigned_msa.split('\n') if line.strip()]
        realigned_msa = "\n".join(realigned_msa[4:])

        return realigned_msa


