import pandas as pd
import Bio.SeqIO as SeqIO

fasta = SeqIO.parse('oxytrijrb310pacbio.faa', 'fasta')
df = pd.DataFrame.from_records([(x.id, str(x.seq)) for x in fasta])
df.columns = ['attr_id', 'sequence']