#!/usr/bin/env python3
"""
Converts a paraest.dat from the old format (all the attractive and
associative parameters written in the file, groups numbered 1..NG) to the
format of GCA v1.14 (PARGR): groups identified by their number in the
database (FOR010-FOR015), plus the parameters given by the user.

    tools/convert_input.py old.dat new.dat [database_dir]

database_dir holds FOR010.DAT-FOR015.DAT (default: the current directory).

Each group of the old file is looked for in FOR010 by (Q, T*) (and g*, g',
g" to break ties). When it is not found, or its database number is already
taken, it gets a database number that the file does not use, and all its
values are written as the user's ("own") values, as one does by hand with a
group that is not in the database. Pure-group values that differ from the
database are written as own values; binary values that differ from what
PARGR would take from the database, too. All the association values are
written as own values (exact, whatever the database has).

The fitted parameters are written with database numbers.

The script reads the file as read_input did (list-directed reads: each READ
starts at a new record and skips blank ones; the rest of the last record is
ignored).
"""
import os
import re
import sys

MISSING = 88.888


# --------------------------------------------------------------------------
# Fortran list-directed reading
# --------------------------------------------------------------------------
class Records:
    def __init__(self, lines):
        self.lines = lines
        self.pos = 0          # next record to read

    def read(self, n, kind=float):
        """READ n values: returns the list of values."""
        vals = []
        while len(vals) < n:
            if self.pos >= len(self.lines):
                raise EOFError('end of file while reading %d values' % n)
            rec = self.pos
            self.pos += 1
            for t in self.lines[rec].replace(',', ' ').split():
                if len(vals) >= n:
                    break
                m = re.fullmatch(r'(\d+)\*(.*)', t)
                if m:
                    vals.extend([conv(m.group(2), kind)] * int(m.group(1)))
                else:
                    vals.append(conv(t, kind))
        return vals[:n]

    def skip_blank(self):
        while self.pos < len(self.lines) and not self.lines[self.pos].split():
            self.pos += 1
        return self.pos

    def read_record(self):
        self.pos += 1
        return self.lines[self.pos - 1]


def conv(t, kind):
    t = t.replace('d', 'e').replace('D', 'E')
    return int(round(float(t))) if kind is int else float(t)


# --------------------------------------------------------------------------
# Database
# --------------------------------------------------------------------------
def read_database(d):
    groups = {}                       # id -> (Q, T*, g*, g', g", name)
    for line in open(os.path.join(d, 'FOR010.DAT'), encoding='latin-1'):
        t = line.split()
        if len(t) >= 6:
            groups[int(t[5])] = (float(t[0]), float(t[1]), float(t[2]), float(t[3]),
                                 float(t[4]), t[6] if len(t) > 6 else 'G%d' % int(t[5]))
    n = len(groups)

    def tokens(name):
        out = []
        for line in open(os.path.join(d, name), encoding='latin-1'):
            out.extend(float(x) for x in line.split())
        return out
    k11, k12, k13 = tokens('FOR011.DAT'), tokens('FOR012.DAT'), tokens('FOR013.DAT')
    kij, kp, alf = {}, {}, {}
    p = 0
    for i in range(1, n + 1):          # row i of FOR011/FOR012: j = i..n
        for j in range(i, n + 1):
            kij[i, j], kp[i, j] = k11[p], k12[p]
            p += 1
    p = 0
    for i in range(1, n + 1):          # row i of FOR013: alpha(i,k), k = 1..n
        for k in range(1, n + 1):
            alf[i, k] = k13[p]
            p += 1
    assoc = {}
    f14 = os.path.join(d, 'FOR014.DAT')
    if os.path.exists(f14):
        for i, line in enumerate(open(f14, encoding='latin-1'), 1):
            t = line.split()
            if len(t) >= 2:
                assoc[i] = (int(t[0]), t[1])
    return groups, kij, kp, alf, assoc


def pargr_binary(db, a, b):
    """k*, k', alpha(a,b), alpha(b,a) that PARGR takes from the database."""
    _, kij, kp, alf, _ = db
    if a == b:
        k, k1 = kij[a, a], 0.0
    else:
        i, j = min(a, b), max(a, b)
        k, k1 = kij[i, j], kp[i, j]
        if abs(k - MISSING) < 1e-3:
            k, k1 = 1.0, 0.0
    al = [0.0 if abs(x - MISSING) < 1e-3 else x for x in (alf[a, b], alf[b, a])]
    return k, k1, al[0], al[1]


def differ(a, b):
    return abs(a - b) > 1e-12 * max(1.0, abs(a), abs(b))


def fmt(x):
    return repr(float(x)) if x != int(x) or abs(x) >= 1e15 else '%.1f' % x


# --------------------------------------------------------------------------
def convert(old, new, dbdir):
    db = read_database(dbdir)
    groups, _, _, _, dbAssoc = db
    raw = open(old, 'rb').read()
    try:
        text, enc = raw.decode('utf-8'), 'utf-8'
    except UnicodeDecodeError:
        text, enc = raw.decode('latin-1'), 'latin-1'
    nl = '\r\n' if '\r\n' in text else '\n'
    lines = text.replace('\r\n', '\n').split('\n')
    if lines and lines[-1] == '':
        lines.pop()
    R = Records(lines)
    warnings = []

    for n in (1, 6, 3, 2, 1, 3):
        R.read(n)
    NC, NG, NGA, NSIG = R.read(4, int)
    ISAT, IGLE, IVE, IUNI, NPRED = R.read(5, int)
    if ISAT:
        idsat = R.read(ISAT, int)
        nsat = R.read(ISAT, int)
        iobj = R.read(ISAT, int)
        for i in range(ISAT):
            for _ in range(nsat[i]):
                R.read(4 if iobj[i] == 81 else 3)

    # ---- old parameter block (replaced)
    blockStart = R.skip_blank()
    Q, Ts = [], []
    for i in range(NG):
        q, t = R.read(2)
        Q.append(q)
        Ts.append(t)
    gs, g1, g2 = R.read(NG), R.read(NG), R.read(NG)
    K = [R.read(NG) for _ in range(NG)]
    K1 = [R.read(NG) for _ in range(NG)]
    A = [R.read(NG) for _ in range(NG)]
    M, E, V = [], {}, {}
    if NGA > 0:
        M = R.read(NGA, int)
        for i in range(NGA):
            for j in range(i, NGA):
                for k in range(M[i]):
                    for l in range(k, M[j]):
                        e, v = R.read(2)
                        for key in ((k, i, l, j), (l, i, k, j), (l, j, k, i), (k, j, l, i)):
                            E[key], V[key] = e, v
    blockEnd = R.pos

    # ---- pure-component properties (kept)
    for i in range(NC):
        R.read(12)
    # ---- group composition (removed: it goes into the new block)
    nuStart = R.skip_blank()
    NU = [R.read(NC, int) for _ in range(NG)]                 # NU[group][comp]
    ONU = [R.read(NC, int) for _ in range(NGA)]
    nuEnd = R.pos

    # ---- data sets (kept)
    for s in range(IGLE):
        ihen = R.read(1, int)[0]
        if ihen == 0:
            ifl = R.read(1, int)[0]
            mez = 2 if ifl < 3 else 3
            nobs = R.read(1 + mez, int)[0]
            nval = 6 if mez == 2 else 10
        elif ihen == 1:
            R.read(1, int)
            nobs = R.read(3, int)[0]
            nval = 5
        elif ihen in (2, 3):
            nobs = R.read(3, int)[0]
            nval = 4 if ihen == 2 else 5
        elif ihen == 5:
            mez = R.read(1, int)[0]
            nobs = R.read(1 + mez, int)[0]
            nval = 7 if mez == 2 else 10
        else:
            raise ValueError('IHEN = %d' % ihen)
        R.read_record()                                        # title
        for _ in range(nobs):
            R.read(nval)

    # ---- fitted parameters (rewritten)
    specStart = R.skip_blank()
    NNP, NNQ, NALFA, NENASS, NRKASS, numTc, NUMDC = R.read(7, int)
    IDG, ITEMP, NIDG, IDIDG = [], [], [], []
    if NNP:
        v = R.read(2 * NNP, int)
        IDG = [(v[2 * i], v[2 * i + 1]) for i in range(NNP)]
        ITEMP = R.read(NNP, int)
        NIDG = R.read(NNP, int)
        for n in NIDG:
            if n:
                v = R.read(2 * n, int)
                IDIDG.append([(v[2 * i], v[2 * i + 1]) for i in range(n)])
            else:
                IDIDG.append([])
    IDQ = R.read(NNQ, int) if NNQ else []
    IDALF, NIDALF, IDIDALF = [], [], []
    if NALFA:
        v = R.read(2 * NALFA, int)
        IDALF = [(v[2 * i], v[2 * i + 1]) for i in range(NALFA)]
        NIDALF = R.read(NALFA, int)
        for n in NIDALF:
            if n:
                v = R.read(2 * n, int)
                IDIDALF.append([(v[2 * i], v[2 * i + 1]) for i in range(n)])
            else:
                IDIDALF.append([])

    def read4(n):
        ids = [tuple(R.read(4, int)) for _ in range(n)]
        nid = R.read(n, int)
        idents = [[tuple(R.read(4, int)) for _ in range(m)] for m in nid]
        return ids, nid, idents
    IDEN, NIDEN, IDIDEN = read4(NENASS) if NENASS else ([], [], [])
    IDRK, NIDRK, IDIDRK = read4(NRKASS) if NRKASS else ([], [], [])
    IDTC = R.read(numTc, int) if numTc else []
    IDDC, ERRDC = (R.read(NUMDC, int), R.read(NUMDC)) if NUMDC else ([], [])
    specEnd = R.pos          # the initial values XX follow (kept)

    # Groups that no component has: PARGR builds the list of groups from the
    # composition lines, so they are listed in component 1 with 0 groups.
    unusedG = [i for i in range(NG) if not any(NU[i])]
    unusedA = [i for i in range(NGA) if not any(ONU[i])]
    if unusedG:
        warnings.append('groups %s are in no component: listed in component 1 with 0 groups'
                        % [i + 1 for i in unusedG])
    if unusedA:
        warnings.append('associating groups %s are in no component: listed in component 1 with 0 groups'
                        % [i + 1 for i in unusedA])

    # ---------------------------------------------------------------------
    # Database numbers of the attractive groups
    # ---------------------------------------------------------------------
    gid = [None] * NG
    used = set()
    for i in range(NG):
        cands = [n for n, g in groups.items()
                 if abs(g[0] - Q[i]) < 1e-4 and abs(g[1] - Ts[i]) < 1e-2]
        cands.sort(key=lambda n: (differ(groups[n][2], gs[i]) + differ(groups[n][3], g1[i])
                                  + differ(groups[n][4], g2[i]), n))
        for n in cands:
            if n not in used:
                gid[i] = n
                used.add(n)
                break
    free = sorted(set(groups) - used, reverse=True)
    unmatched = {i + 1 for i in range(NG) if gid[i] is None}
    for i in range(NG):
        if gid[i] is None:
            gid[i] = free.pop(0)
            used.add(gid[i])
            warnings.append('group %d (Q = %g, T* = %g) is not in the database: it takes number %d '
                            'with own values' % (i + 1, Q[i], Ts[i], gid[i]))
    name = {}
    for i in range(NG):
        name[i] = groups[gid[i]][5][:10] if (i + 1) not in unmatched else 'G%d' % (i + 1)

    # Associating groups: a free number among those with the same sites
    aid = [None] * NGA
    usedA = set()
    for i in range(NGA):
        cands = [n for n, (s, _) in sorted(dbAssoc.items()) if s == M[i] and n not in usedA]
        if not cands:
            cands = [n for n in sorted(dbAssoc) if n not in usedA]
        aid[i] = cands[0]
        usedA.add(aid[i])
    if NGA:
        warnings.append('associating groups: database numbers %s chosen only by their number of '
                        'sites; all their values are own values' % aid)

    # ---------------------------------------------------------------------
    # New block
    # ---------------------------------------------------------------------
    out = []
    pure = []
    for i in range(NG):
        g = groups[gid[i]]
        vals = (Q[i], Ts[i], gs[i], g1[i], g2[i])
        if any(differ(a, b) for a, b in zip(vals, g[:5])):
            pure.append('%d  %s  %s  %s  %s  %s  %s' % ((gid[i],) + tuple(fmt(x) for x in vals) + (name[i],)))
    binary = []
    for i in range(NG):
        for j in range(i, NG):
            if differ(K[i][j], K[j][i]) or differ(K1[i][j], K1[j][i]):
                warnings.append('k*/k\' of groups %d-%d are not symmetric: (%d,%d) is used'
                                % (i + 1, j + 1, i + 1, j + 1))
            want = (K[i][j], K1[i][j], A[i][j], A[j][i])
            have = pargr_binary(db, gid[i], gid[j])
            if any(differ(a, b) for a, b in zip(want, have)):
                binary.append('%d  %d  %s  %s  %s  %s' % ((gid[i], gid[j]) + tuple(fmt(x) for x in want)))
    out.append('%d  %d    ! NownII, NownIJ (attractive groups)' % (len(pure), len(binary)))
    for c in range(NC):
        # in increasing order of the database number: group_reading (PARGR)
        # needs it
        pairs = sorted((gid[i], NU[i][c]) for i in range(NG) if NU[i][c] != 0 or (c == 0 and i in unusedG))
        toks = ['%d %d' % p for p in pairs]
        pad = 2 * NG - 2 * len(pairs)
        out.append('  '.join(toks) + ('  %d*0' % pad if pad else '') + '    ! component %d' % (c + 1))
    out.extend(pure)
    out.extend(binary)
    if NGA > 0:
        sites = ['%d  %d  A%d' % (aid[i], M[i], i + 1) for i in range(NGA)]
        pairs = []
        for i in range(NGA):
            for j in range(i, NGA):
                for k in range(M[i]):
                    for l in range(M[j]):
                        if i == j and l < k:
                            continue
                        # Value used by the old program (fakePARGR): enAss(l,j,k,i)
                        # for site k of group i and site l of group j (i <= j).
                        # PARGR wants the group with the lower number first.
                        kk, ii, ll, jj = (k, i, l, j) if aid[i] <= aid[j] else (l, j, k, i)
                        pairs.append('%d  %d  %d  %d  %s  %s' % (kk + 1, aid[ii], ll + 1, aid[jj],
                                                              fmt(E.get((l, j, k, i), 0.0)),
                                                              fmt(V.get((l, j, k, i), 0.0))))
        out.append('%d  %d    ! NownII, NownIJ (associating groups)' % (len(sites), len(pairs)))
        for c in range(NC):
            pr = sorted((aid[i], ONU[i][c]) for i in range(NGA) if ONU[i][c] != 0 or (c == 0 and i in unusedA))
            toks = ['%d %d' % p for p in pr]
            pad = 2 * NGA - 2 * len(pr)
            out.append('  '.join(toks) + ('  %d*0' % pad if pad else '') + '    ! component %d' % (c + 1))
        out.extend(sites)
        out.extend(pairs)

    # ---------------------------------------------------------------------
    # New fitted-parameter specification (database numbers)
    # ---------------------------------------------------------------------
    G = lambda i: gid[i - 1]
    GA = lambda i: aid[i - 1]
    spec = ['%d %d %d %d %d %d %d    ! NNP NNQ NALFA NENASS NRKASS numTc NUMDC'
            % (NNP, NNQ, NALFA, NENASS, NRKASS, numTc, NUMDC)]
    if NNP:
        spec.append('  '.join('%d %d' % (G(a), G(b)) for a, b in IDG) + '    ! g/k (database numbers)')
        spec.append(' '.join(str(x) for x in ITEMP))
        spec.append(' '.join(str(x) for x in NIDG))
        for lst in IDIDG:
            if lst:
                spec.append('  '.join('%d %d' % (G(a), G(b)) for a, b in lst))
    if NNQ:
        spec.append(' '.join(str(G(i)) for i in IDQ) + '    ! areas')
    if NALFA:
        spec.append('  '.join('%d %d' % (G(a), G(b)) for a, b in IDALF) + '    ! alpha')
        spec.append(' '.join(str(x) for x in NIDALF))
        for lst in IDIDALF:
            if lst:
                spec.append('  '.join('%d %d' % (G(a), G(b)) for a, b in lst))
    for ids, nid, idents, what in ((IDEN, NIDEN, IDIDEN, 'association energy'),
                                   (IDRK, NIDRK, IDIDRK, 'association volume')):
        if not ids:
            continue
        for k, i, l, j in ids:
            spec.append('%d %d %d %d    ! %s' % (k, GA(i), l, GA(j), what))
        spec.append(' '.join(str(x) for x in nid))
        for lst in idents:
            for k, i, l, j in lst:
                spec.append('%d %d %d %d' % (k, GA(i), l, GA(j)))
    if IDTC:
        spec.append(' '.join(str(x) for x in IDTC))
    if IDDC:
        spec.append(' '.join(str(x) for x in IDDC))
        spec.append(' '.join(fmt(x) for x in ERRDC))

    result = (lines[:blockStart] + out + [''] + lines[blockEnd:nuStart] + lines[nuEnd:specStart]
              + spec + lines[specEnd:])
    open(new, 'w', encoding=enc, newline='').write(nl.join(result) + nl)
    return warnings, gid, aid


if __name__ == '__main__':
    w, gid, aid = convert(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else '.')
    print('attractive groups ->', gid)
    if aid:
        print('associating groups ->', aid)
    for x in w:
        print('WARNING:', x)
