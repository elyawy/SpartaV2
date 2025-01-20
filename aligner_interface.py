import os, pathlib, tempfile, functools
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import ClustalwCommandline


aligner_params_dict = {
    "MAFFT": functools.partial(MafftCommandline,thread=16 , globalpair=True, maxiterate=1000),
    "PRANK": functools.partial(PrankCommandline, f="fasta", once=False),
    "CLUSTALW": functools.partial(ClustalwCommandline, output="FASTA")
}


class Aligner:
    def __init__(self, aligner: str) -> None:
        self._aligner_name = aligner
        self._aligner = aligner_params_dict[aligner]        
    
    def set_input_file(self, file_path, tree_file=None):
        self._input_file = str(file_path)
        self._output_name = str(file_path).split(".")[0]
        if self._aligner_name == "MAFFT":
            self._aligner_cmd = self._aligner(input=self._input_file)
        if self._aligner_name == "PRANK":
            self._aligner_cmd = self._aligner(d=self._input_file, o=self._output_name, t=str(tree_file))
        if self._aligner_name == "CLUSTALW":
            self._aligner_cmd = self._aligner(infile=self._input_file, outfile=f"{self._output_name}_realigned.fasta")


    def get_realigned_msa(self) -> str:
        # print(self._aligner_cmd)
        if self._aligner_name == "MAFFT":
            realigned_msa, stderr = self._aligner_cmd()
        if self._aligner_name == "PRANK":
            realigned_msa, stderr = self._aligner_cmd()
            with open(f"{self._output_name}.best.fas",'r') as f:
                realigned_msa = f.read()
        if self._aligner_name == "CLUSTALW":
            realigned_msa, stderr = self._aligner_cmd()
            with open(f"{self._output_name}_realigned.fasta",'r') as f:
                realigned_msa = f.read()


        return realigned_msa


