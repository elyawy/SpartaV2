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
            self._aligner_cmd = ["mafft", "--globalpair", "--maxiterate", "0", self._input_file]


    def get_realigned_msa(self) -> str:
        # print(self._aligner_cmd)
        if self._aligner_name == "MAFFT":
            result = subprocess.run(self._aligner_cmd, capture_output=True, text=True)
            realigned_msa, stderr = result.stdout, result.stderr

        return realigned_msa


