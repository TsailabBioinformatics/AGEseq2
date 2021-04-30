"""
Utility class for AgeseqIO to handle different types of user files.

File errors  propagate as ValueErrors up to calling functions.
"""

from pathlib import Path, PurePath
from Bio import SeqIO
import gzip
import os


class ReadsContextManager:
    """ Custom ContextManager class to allow use with gzip versus regular files.  If uncompressed, this behaves
    exactly like a normal with open(file, mode).  If the file is compressed, it will use gzip to open the file and
    return the context. """
    def __init__(self, filename, mode, compressed=False):
        self.filename = filename
        self.mode = mode
        self.compressed = compressed
        self.file = None

    def __enter__(self):
        if self.compressed:
            self.file = gzip.open(self.filename, self.mode)
            return self.file
        else:
            self.file = open(self.filename, self.mode)
            return self.file

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.file.close()


class ReadsHandler:

    def __init__(self, working_dir=Path.cwd()):
        self.working_dir = working_dir

    @property
    def reads_dir(self):
        """This property is set outside of init to ensure path and OS consistency.
        It also allows test files without crowding or using separate unit tests
        outside this module.
        """
        # This enforces path consistency in case a user manually inputs a String-like path.
        if not isinstance(self.working_dir, PurePath):
            self.working_dir = Path(self.working_dir)

        # Running the code form main injects a RH_TEST=1 environmental variable which is caught here to alter the
        # path to the test files.
        if os.environ.get('RH_TEST'):
            return self.working_dir / 'test_files'

        return self.working_dir / 'reads'

    @staticmethod
    def _is_compressed(file_loc):
        """ Returns whether the file is actually gzipped or not. """
        with open(file_loc, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'  # magic number

    @staticmethod
    def _get_file_type(file_id):
        """ Biopython needs the filetype as fasta or fastq. This attempts to
        extract it.
        """
        if 'fastq' in file_id:
            file_type = 'fastq'
        elif 'fasta' in file_id:
            file_type = 'fasta'
        else:
            raise ValueError("Please ensure file extension ends with .fasta .fastq, .fasta.gz, or .fastq.gz")

        return file_type

    def handle(self, file_id):
        """ Receives a file_id and performs fasta/fastq and gzip checks.  If
        everything is fine, this uses BioPython to parse the file and return
        a dictionary of the file's contents.
        """
        file_loc = self.reads_dir / file_id

        if not file_loc.is_file():
            raise ValueError("Error finding file {}".format(file_loc))

        mode = 'r'
        compressed = False

        # Check compression status to not waste time or error on falsely compressed
        if self._is_compressed(file_loc):
            mode = 'rt'
            compressed = True

        # Bio can silently return an empty dictionary so at every step it is important to bind {} to None for
        # consistency.
        with ReadsContextManager(file_loc, mode, compressed) as read_file:
            file_type = self._get_file_type(file_id)
            try:
                seq_parser = SeqIO.to_dict(SeqIO.parse(read_file, file_type))
            except ValueError as ex:
                # Propagate error up to calling function
                raise

        if seq_parser == {}:
            # return None
            raise ValueError("Error parsing!  BioPython returned empty dict.")

        return seq_parser


if __name__ == "__main__":
    import unittest
    os.environ['RH_TEST'] = '1'
    rh = ReadsHandler(os.getcwd())
    contents = rh.handle('SBHP73_S73_L001_R1_001.fastq')
    # print(contents)

    class ReadsHandlerTests(unittest.TestCase):
        @staticmethod
        def test_good_file():
            """ Ensures a proper fastq file has return information. """
            assert rh.handle('SBHP72_S72_L001_R2_001.fastq') is not None

        @staticmethod
        def test_missing_compression_label():
            """ Checks a gzipped fasta file is properly decompressed. """
            assert rh.handle('FASTA_MISSING_GZ.fasta') is not None

        @staticmethod
        def test_false_compression_label():
            """ Checks a false compression label is not decompressed. """
            assert rh.handle('FASTA_FALSE_GZ.fasta.gz') is not None

        def test_bad_filename(self):
            """ Checks an improper file throws a ValueError. """
            self.assertRaises(
                ValueError,
                rh.handle,
                'badfile.txt'
            )

        def test_bad_input(self):
            """ Ensures errors from BioPython are propagated up to callers. """
            self.assertRaises(
                ValueError,
                rh.handle,
                'badfile.fastq'
            )


    unittest.main()