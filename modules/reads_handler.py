from pathlib import Path, PurePath
from Bio import SeqIO
import gzip
import os


class ReadsHandler:

	def __init__(self, working_dir=Path.cwd()):
		self.working_dir = working_dir

	@property
	def reads_dir(self):
		'''This property is set outside of init to ensure path and OS consistency.
		It also allows test files without crowding or using separate unit tests
		outside this module.
		'''
		if not isinstance(self.working_dir, PurePath):
			self.working_dir = Path(self.working_dir)

		if os.environ.get('RH_TEST'):
			return self.working_dir / 'test_files'

		return self.working_dir / 'reads'

	def _is_compressed(self, file_loc):
		''' Returns whether the file is actually gzipped or not. '''
		with open(file_loc, 'rb') as f:
			return f.read(2) == b'\x1f\x8b'  # magic number

	def _get_file_type(self, file_id):
		''' Biopython needs the filetype as fasta or fastq. This attempts to 
		extract it.
		'''
		if 'fastq' in file_id:
			file_type = 'fastq'
		elif 'fasta' in file_id:
			file_type = 'fasta'
		else:
			raise ValueError("Please ensure file extension ends with .fasta .fastq, .fasta.gz, or .fastq.gz")

		return file_type

	def handle(self, file_id):
		''' Receives a file_id and performs fasta/fastq and gzip checks.  If
		everything is fine, this uses BioPython to parse the file and return 
		a dictionary of the file's contents.
		'''
		file_loc = self.reads_dir / file_id

		if not file_loc.is_file():
			raise ValueError("Error finding file {}".format(file_loc))

		# # Check compression status to not waste time or error on falsely compressed
		if self._is_compressed(file_loc):
			file_loc = gzip.open(file_loc, 'rt')
		else:
			file_loc = open(file_loc, 'r')
		
		# Try and parse the file - this will error out if not named correctly.
		# If the file is named correctly and opened but incorrect, BioPython
		# silently returns an empty dictionary
		file_type = self._get_file_type(file_id)
		try:
			seq_parser = SeqIO.to_dict(SeqIO.parse(file_loc, file_type))
		except ValueError as e:
			print("!! ERROR PARSING FILE.  BioPython rejected contents. !!")
			return None
		finally:
			file_loc.close()
		
		# Bind a null return since BioPython can pass an empty dictionary
		if seq_parser == {}:
			file_loc.close()
			raise ValueError("Error parsing!  BioPython returned empty dict.")
		
		file_loc.close()
		return seq_parser


if __name__ == "__main__":
	import unittest
	os.environ['RH_TEST'] = '1'

	rh = ReadsHandler(os.getcwd())
	
	class ReadsHandlerTests(unittest.TestCase):
		def test_good_file(self):
			assert rh.handle('SBHP72_S72_L001_R2_001.fastq') != None
		
		def test_missing_compression_label(self):
			assert rh.handle('FASTA_MISSING_GZ.fasta') != None
				
		def test_false_compression_label(self):
			assert rh.handle('FASTA_FALSE_GZ.fasta.gz') != None
		
		def test_bad_filename(self):
			self.assertRaises(
				ValueError,
				rh.handle,
				'badfile.txt'
			)
		
		def test_bad_input(self):
			assert rh.handle('badfile.fastq') == None
			
	
	unittest.main()
	