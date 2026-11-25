#!/usr/bin/env python3
"""
Converts a paraest.dat from the old format (all the attractive and
associative parameters written in the file, groups numbered 1..NG) to the
format of GCA v1.14 (PARGR): groups identified by their number in the
database (FOR010-FOR015), plus the parameters given by the user.

    tools/convert_input.py old.dat new.dat [database_dir]
    tools/convert_input.py --components v114.dat new.dat

database_dir holds FOR010.DAT-FOR015.DAT (default: the current directory).

Second stage (applied always; alone with --components, for files already
in the v1.14 format): the line of each component "Tc Pc Tsat Psat dc Tfu
DHfu DCpfu TtrSS DHss DVfu Pfu" is replaced by the input of PARMOL (unit
option 1 = atm; per component its name, "Tc Pc omega Tsat Psat dc" with
omega = 0, "M c A B C D E F G" all 0) followed by the melting block, one
line per component "Tfus Pfus DHfus DVfus ACpfus BCpfus NTSS" (BCpfus = 0;
NTSS = 1 when the old line had a solid-solid transition, followed by its
line "Ttr DHtr"). Components are named comp1, comp2, ... (the comment of
the old line, if any, goes to the line of the name).

Each group of the old file is looked for in FOR010 by (Q, T*) (and g*, g',
g" to break ties). When it is not found, or its database number is already
taken, it gets a database number that the file does not use, and all its
values are written as the user's ("own") values, as one does by hand with a
group that is not in the database. Pure-group values that differ from the
database are written as own values; binary values that differ from what
PARGR would take from the database, too; and the association energies and
volumes that differ from FOR015 (a pair missing from the database is 0).

Each associating group is identified by the name in the comment of its
composition line ("!asoc AOH", against the short name or the description of
FOR014), else by its association energies and volumes (with itself and
with the groups already identified) against FOR015; a group found neither
way takes a free number, with its sites, name and parameters written as own
values.

The fitted parameters are written with the position of each group in the
new system: groups in increasing order of their database number (column Pos
of the parameter tables in paraest.out).

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
                assoc[i] = (int(t[0]), t[1], ' '.join(t[3:]))   # sites, short name, description
    # FOR015, in the order PARGR reads it: for each group a and site k, the
    # groups b >= a and their sites l (l >= k when b = a). 88.888 = missing,
    # taken as 0 (as PARGR does).
    assocPair = {}
    f15 = os.path.join(d, 'FOR015.DAT')
    if os.path.exists(f15):
        recs = [l.split() for l in open(f15, encoding='latin-1') if l.split()]
        p = 0
        n = len(assoc)
        for a in range(1, n + 1):
            for k in range(1, assoc[a][0] + 1):
                for b in range(a, n + 1):
                    for l in range(k if b == a else 1, assoc[b][0] + 1):
                        e, v = float(recs[p][0]), float(recs[p][1])
                        p += 1
                        if abs(e - MISSING) < 1e-3:
                            e, v = 0.0, 0.0
                        assocPair[k, a, l, b] = (e, v)
                        assocPair[l, b, k, a] = (e, v)
    return groups, kij, kp, alf, assoc, assocPair


def pargr_binary(db, a, b):
    """k*, k', alpha(a,b), alpha(b,a) that PARGR takes from the database."""
    _, kij, kp, alf, _, _ = db
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


def norm(name):
    """Name without case, '-', '_', spaces...: 'o-AO-' -> 'oao'."""
    return re.sub(r'[^0-9a-z]', '', name.lower())


ALIASES = {'ar': 'aring'}


def assoc_ids(NGA, M, E, V, onuName, dbAssoc, dbPair):
    """Database numbers of the associating groups of the old file.

    1. By name: the comment of the composition line ("!asoc AOH") against the
       short name or the description of FOR014.
    2. By value: the association energies and volumes of the group with
       itself and with the groups already identified, against FOR015;
       repeated while it identifies groups (the unique candidate with the
       most coincident nonzero pairs).
    3. The rest take the highest free numbers, with own values.
    Returns the numbers and how each one was found ('name', 'value', 'own').
    """
    aid, how = [None] * NGA, [None] * NGA
    names = {}
    for n, (_, short, desc) in dbAssoc.items():
        for x in (short, desc):
            names.setdefault(norm(x), set()).add(n)
    used = set()
    for i in range(NGA):
        key = ALIASES.get(norm(onuName[i]), norm(onuName[i]))
        cands = names.get(key, set()) - used if key else set()
        if len(cands) == 1:
            aid[i], how[i] = cands.pop(), 'name'
            used.add(aid[i])

    def close(a, b):
        return abs(a - b) <= 1e-6 * max(1.0, abs(a), abs(b))

    def score(i, n):
        hits = 0
        for j in range(NGA):
            m = n if j == i else aid[j]
            if m is None:
                continue
            for k in range(M[i]):
                for l in range(M[j]):
                    de, dv = dbPair.get((k + 1, n, l + 1, m), (0.0, 0.0))
                    if (de or dv) and close(E.get((k, i, l, j), 0.0), de) and close(V.get((k, i, l, j), 0.0), dv):
                        hits += 1
        return hits
    while True:
        best = []
        for i in range(NGA):
            if aid[i] is not None:
                continue
            for n in dbAssoc:
                if n not in used and dbAssoc[n][0] == M[i]:
                    h = score(i, n)
                    if h:
                        best.append((h, i, n))
        if not best:
            break
        h, i, n = max(best)
        if sum(1 for b in best if b[0] == h and (b[1] == i or b[2] == n)) > 1:
            break                              # ambiguous: left for step 3
        aid[i], how[i] = n, 'value'
        used.add(n)
    free = sorted(set(dbAssoc) - used, reverse=True)
    for i in range(NGA):
        if aid[i] is None:
            aid[i], how[i] = free.pop(0), 'own'
    return aid, how


# --------------------------------------------------------------------------
def convert(old, new, dbdir):
    db = read_database(dbdir)
    groups, _, _, _, dbAssoc, dbPair = db
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
    def comment_name():
        # the name of the group, if the line ends with a comment ("!asoc AOH")
        c = lines[R.pos - 1].split('!', 1)
        name = c[1].split()[-1] if len(c) > 1 and c[1].split() else ''
        return name if name.isascii() and norm(name) else ''
    NU, nuName = [], []                                       # NU[group][comp]
    for _ in range(NG):
        NU.append(R.read(NC, int))
        nuName.append(comment_name())
    ONU, onuName = [], []
    for _ in range(NGA):
        ONU.append(R.read(NC, int))
        onuName.append(comment_name())
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
        name[i] = groups[gid[i]][5][:10] if (i + 1) not in unmatched else (nuName[i] or 'G%d' % (i + 1))[:10]

    aid, how = assoc_ids(NGA, M, E, V, onuName, dbAssoc, dbPair)
    for i in range(NGA):
        if how[i] == 'own':
            warnings.append('associating group %d (%s, %d sites) is not in the database: it takes number %d '
                            'with own sites and association parameters' % (i + 1, onuName[i] or '?', M[i], aid[i]))
        elif how[i] == 'name' and dbAssoc[aid[i]][0] != M[i]:
            warnings.append('associating group %d (%s): %d sites, %d in the database (own sites written)'
                            % (i + 1, onuName[i], M[i], dbAssoc[aid[i]][0]))

    # ---------------------------------------------------------------------
    # New block
    # ---------------------------------------------------------------------
    out = []
    pure = []
    for i in range(NG):
        g = groups[gid[i]]
        vals = (Q[i], Ts[i], gs[i], g1[i], g2[i])
        if any(differ(a, b) for a, b in zip(vals, g[:5])):
            pure.append('%d  %s  %s  %s  %s  %s  %s' % ((gid[i],) + tuple(fmt(x) for x in vals) + (name[i],))
                        + ('    ! not in the database (old group %d)' % (i + 1) if (i + 1) in unmatched else ''))
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
        # own sites only when the database has another number of sites
        sites = ['%d  %d  %s' % (aid[i], M[i], (onuName[i] or 'A%d' % (i + 1))[:10])
                 + ('    ! not in the database (old associating group %d)' % (i + 1) if how[i] == 'own' else '')
                 for i in range(NGA) if how[i] == 'own' or dbAssoc.get(aid[i], (None,))[0] != M[i]]
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
                        e, v = E.get((l, j, k, i), 0.0), V.get((l, j, k, i), 0.0)
                        # only what differs from the database (missing = 0)
                        de, dv = dbPair.get((kk + 1, aid[ii], ll + 1, aid[jj]), (0.0, 0.0))
                        if differ(e, de) or differ(v, dv):
                            pairs.append('%d  %d  %d  %d  %s  %s' % (kk + 1, aid[ii], ll + 1, aid[jj],
                                                                  fmt(e), fmt(v)))
        out.append('%d  %d    ! NownII, NownIJ (associating groups)' % (len(sites), len(pairs)))
        for c in range(NC):
            pr = sorted((aid[i], ONU[i][c]) for i in range(NGA) if ONU[i][c] != 0 or (c == 0 and i in unusedA))
            toks = ['%d %d' % p for p in pr]
            pad = 2 * NGA - 2 * len(pr)
            out.append('  '.join(toks) + ('  %d*0' % pad if pad else '') + '    ! component %d' % (c + 1))
        out.extend(sites)
        out.extend(pairs)

    # ---------------------------------------------------------------------
    # New fitted-parameter specification (positions in the new system)
    # ---------------------------------------------------------------------
    # position in the new system: PARGR sorts the groups by database number
    G = lambda i: sorted(gid).index(gid[i - 1]) + 1
    GA = lambda i: sorted(aid).index(aid[i - 1]) + 1
    spec = ['%d %d %d %d %d %d %d    ! NNP NNQ NALFA NENASS NRKASS numTc NUMDC'
            % (NNP, NNQ, NALFA, NENASS, NRKASS, numTc, NUMDC)]
    if NNP:
        spec.append('  '.join('%d %d' % (G(a), G(b)) for a, b in IDG) + '    ! g/k (group positions)')
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
    result = convert_components(result)
    open(new, 'w', encoding=enc, newline='').write(nl.join(result) + nl)
    return warnings, gid, aid


def convert_components(lines):
    """v1.14 format -> PARMOL input + melting block (see the docstring)."""
    R = Records(lines)
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
    # PARGR block: attractive part, then (NGA > 0) associating part
    nII, nIJ = R.read(2, int)
    for _ in range(NC):
        R.read(2 * NG)
    for _ in range(nII):
        R.skip_blank(); R.read_record()          # id, 5 values, name
    for _ in range(nIJ):
        R.read(6)
    if NGA > 0:
        nII, nIJ = R.read(2, int)
        for _ in range(NC):
            R.read(2 * NGA)
        for _ in range(nII):
            R.skip_blank(); R.read_record()      # id, sites, name
        for _ in range(nIJ):
            R.read(6)
    start = R.skip_blank()
    comp, note = [], []
    for _ in range(NC):
        comp.append(R.read(12))
        c = lines[R.pos - 1].split('!', 1)
        note.append('    !' + c[1].rstrip() if len(c) > 1 and c[1].strip() else '')
    end = R.pos
    out = ['1    ! PARMOL: units (1: atm)']
    for i, (Tc, Pc, Tsat, Psat, dc, Tfu, DHfu, DCpfu, Ttr, DHss, DVfu, Pfu) in enumerate(comp):
        out.append('comp%d' % (i + 1) + note[i])
        out.append('  '.join(fmt(x) for x in (Tc, Pc, 0.0, Tsat, Psat, dc)) + '    ! Tc Pc omega Tsat Psat dc')
        out.append('0 0 0 0 0 0 0 0 0    ! M c A B C D E F G')
    out.append('')
    for i, (Tc, Pc, Tsat, Psat, dc, Tfu, DHfu, DCpfu, Ttr, DHss, DVfu, Pfu) in enumerate(comp):
        ntss = 1 if Ttr != 0 or DHss != 0 else 0
        out.append('  '.join(fmt(x) for x in (Tfu, Pfu, DHfu, DVfu, DCpfu, 0.0)) + '  %d' % ntss
                   + '    ! Tfus Pfus DHfus DVfus ACpfus BCpfus NTSS (comp%d)' % (i + 1))
        if ntss:
            out.append('  '.join(fmt(x) for x in (Ttr, DHss)) + '    ! Ttr DHtr')
    return lines[:start] + out + lines[end:]


if __name__ == '__main__' and sys.argv[1] == '--components':
    raw = open(sys.argv[2], 'rb').read()
    try:
        text, enc = raw.decode('utf-8'), 'utf-8'
    except UnicodeDecodeError:
        text, enc = raw.decode('latin-1'), 'latin-1'
    nl = '\r\n' if '\r\n' in text else '\n'
    lines = text.replace('\r\n', '\n').split('\n')
    if lines and lines[-1] == '':
        lines.pop()
    open(sys.argv[3], 'w', encoding=enc, newline='').write(nl.join(convert_components(lines)) + nl)
elif __name__ == '__main__':
    w, gid, aid = convert(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else '.')
    print('attractive groups ->', gid)
    if aid:
        print('associating groups ->', aid)
    for x in w:
        print('WARNING:', x)
