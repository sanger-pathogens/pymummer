class Error(Exception):
    pass


class Snp:
    def __init__(self, line):
        # Without the -C option to show-snps, looks like this:
        # [P1] [SUB] [SUB]  [P2] [BUFF] [DIST] [LEN R] [LEN Q] [FRM]   [TAGS]
        # 187  A     C      269  187    187    654     853     1       1   ref_name  qry_name

        # With the -C option to show-snps, looks like this:
        # [P1] [SUB] [SUB] [P2] [BUFF] [DIST] [R] [Q] [LEN R] [LEN Q] [FRM] [TAGS]
        # 187  A     C     269  187    187    0   0   654     853     1     1   ref_name  qry_name
        try:
            l = line.rstrip().split("\t")
            self.ref_pos = int(l[0]) - 1
            self.ref_base = l[1]
            self.qry_base = l[2]
            self.qry_pos = int(l[3]) - 1
            self.ref_length = int(l[-6])
            self.qry_length = int(l[-5])
            self.reverse = {"1": False, "-1": True}[l[-3]]
            self.ref_name = l[-2]
            self.qry_name = l[-1]
        except:
            raise Error(
                "Error constructing pymummer.snp.Snp from mummer show-snps output at this line:\n"
                + line
            )

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __str__(self):
        return "\t".join(
            [
                str(x)
                for x in [
                    self.ref_pos + 1,
                    self.ref_base,
                    self.qry_base,
                    self.qry_pos + 1,
                    self.ref_length,
                    self.qry_length,
                    "-1" if self.reverse else "1",
                    self.ref_name,
                    self.qry_name,
                ]
            ]
        )
