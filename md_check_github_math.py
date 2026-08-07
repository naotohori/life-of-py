#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Check LaTeX math in Markdown for constructs GitHub renders incorrectly.

GitHub runs Markdown processing over a file *before* handing the math to
KaTeX. That silently rewrites some LaTeX, so an expression that is perfectly
valid on its own can render as garbage, or not at all, once it is in a .md
file. This tool finds those cases.

Three failure classes are detected statically (no dependencies, no network):

  escape    A backslash before ASCII punctuation is a Markdown escape and the
            backslash is removed. So \\; \\, \\! reach KaTeX as bare ; , !
            and are typeset as literal punctuation:
                E \\;=\\; U      renders as    E ;=; U
            Fix: delete them, or use an alphabetic macro such as \\thinspace,
            \\quad or \\qquad, which Markdown does not touch.

  textmode  A consequence of the above: \\texttt{ef\\_tmp} loses its backslash,
            leaving a bare _ inside \\texttt{}, which is text mode. KaTeX then
            fails outright with "'_' allowed only in math mode".
            Fix: \\texttt{ef}_\\texttt{tmp}, a normal math-mode subscript.

  emphasis  An underscore preceded by punctuation (typically } or ) ) can open
            Markdown emphasis. Markdown pairs it with a later underscore,
            consumes both, and wraps the span in <em>, wrecking the math and
            usually the rest of the paragraph:
                $\\mathbf{R}_{l,t}$ ... $U_\\alpha$   ->   <em> ate both _
            Fix: keep every underscore after an alphanumeric, e.g. \\mathbf R
            rather than \\mathbf{R}, and \\rho_i^{(b)} rather than \\rho^{(b)}_i.
            Display ($$) blocks are immune; only inline math is at risk.

Two optional deeper checks:

  --katex   Parse every expression with KaTeX itself, after emulating the
            Markdown unescaping, catching any other LaTeX error. Needs node
            with the katex module resolvable (npm install -g katex).
  --api     Ground truth: render via GitHub's own Markdown API and report any
            math span that came back split by <em>/<strong>. Needs network;
            unauthenticated use is rate limited (~60 requests/hour).

Exit status is 1 if anything was reported, so it works in a pre-commit hook.
"""

import argparse
import html
import json
import os
import re
import subprocess
import sys
import tempfile
import urllib.request

ASCII_PUNCT = r'!"#$%&\'()*+,\-./:;<=>?@\[\\\]^_`{|}~'
ESCAPE_RE = re.compile(r'\\([' + ASCII_PUNCT + r'])')
TEXT_CMDS = ('text', 'texttt', 'textrm', 'textbf', 'textit', 'textsf', 'mbox', 'hbox')


class Finding:
    def __init__(self, path, line, kind, detail, snippet, fix):
        self.path, self.line, self.kind = path, line, kind
        self.detail, self.snippet, self.fix = detail, snippet, fix

    def report(self, colour):
        red, dim, off = ('\033[31m', '\033[2m', '\033[0m') if colour else ('', '', '')
        print(f'{red}{self.path}:{self.line}{off}  [{self.kind}] {self.detail}')
        print(f'    {dim}{self.snippet}{off}')
        print(f'    fix: {self.fix}')


def strip_code(text):
    """Blank out fenced blocks and inline code spans, preserving line numbers."""
    def blank(m):
        return re.sub(r'[^\n]', ' ', m.group(0))
    text = re.sub(r'^([ \t]*)(```|~~~).*?^\1\2[^\n]*$', blank, text, flags=re.S | re.M)
    text = re.sub(r'(?<!`)(`+)(?!`)(.+?)(?<!`)\1(?!`)', blank, text, flags=re.S)
    return text


def find_math(text):
    """Yield (start, end, tex, display) for each math span. Input is code-stripped."""
    spans = []
    for m in re.finditer(r'^[ \t]*\$\$[ \t]*\n(.*?)\n[ \t]*\$\$[ \t]*$', text, re.S | re.M):
        spans.append((m.start(), m.end(), m.group(1), True))
    masked = list(text)
    for s, e, _, _ in spans:
        for i in range(s, e):
            masked[i] = ' '
    masked = ''.join(masked)
    for m in re.finditer(r'(?<!\$)\$(?!\s)([^$\n]+?)(?<!\s)\$(?!\$)', masked):
        spans.append((m.start(), m.end(), m.group(1), False))
    return sorted(spans)


def unescape(tex):
    """What Markdown hands to KaTeX after its own escape processing."""
    return ESCAPE_RE.sub(r'\1', tex)


def can_open_emphasis(s, i):
    """Can the underscore at s[i] open Markdown emphasis? (left-flanking rule)"""
    nxt = s[i + 1] if i + 1 < len(s) else ' '
    prv = s[i - 1] if i > 0 else ' '
    if nxt.isspace():
        return False
    if not nxt.isalnum():            # followed by punctuation: needs space/punct before
        return not prv.isalnum()
    return not prv.isalnum()         # '_' between alphanumerics never opens


def line_of(text, pos):
    return text.count('\n', 0, pos) + 1


def check_file(path, want_katex, want_api):
    raw = open(path, encoding='utf-8').read()
    text = strip_code(raw)
    findings, exprs = [], []

    for start, _end, tex, display in find_math(text):
        line = line_of(text, start)
        exprs.append((line, tex, display))
        seen = unescape(tex)

        for m in ESCAPE_RE.finditer(tex):
            seq = m.group(0)
            findings.append(Finding(
                path, line, 'escape',
                f'{seq!r} loses its backslash and reaches KaTeX as {m.group(1)!r}',
                tex.strip()[:100],
                'drop it, or use \\thinspace / \\quad / \\qquad instead'))
            break

        for cmd in TEXT_CMDS:
            for m in re.finditer(r'\\' + cmd + r'\{([^{}]*)\}', seen):
                if '_' in m.group(1) or '^' in m.group(1):
                    findings.append(Finding(
                        path, line, 'textmode',
                        f'"_" or "^" inside \\{cmd}{{}} (text mode); KaTeX rejects it',
                        m.group(0), r'use a math-mode subscript, e.g. \texttt{a}_\texttt{b}'))

        if not display:
            for i, ch in enumerate(tex):
                if ch == '_' and can_open_emphasis(tex, i):
                    findings.append(Finding(
                        path, line, 'emphasis',
                        f'"_" after {tex[i-1]!r} can open Markdown emphasis and be eaten',
                        tex.strip()[:100],
                        'keep every "_" after a letter or digit '
                        r'(\mathbf R not \mathbf{R}; \rho_i^{(b)} not \rho^{(b)}_i)'))
                    break

    if want_katex:
        findings += katex_check(path, exprs)
    if want_api:
        findings += api_check(path, raw)
    return findings, len(exprs)


NODE_SCRIPT = r'''
const katex = require('katex');
let out = [];
for (const it of JSON.parse(process.argv[1])) {
  try { katex.renderToString(it.tex, {displayMode: it.display, throwOnError: true}); }
  catch (e) { out.push({line: it.line, msg: e.message.split('\n')[0]}); }
}
console.log(JSON.stringify(out));
'''


VENDOR = os.path.expanduser('~/.local/share/katex/node_modules')
_katex_env = None   # cached: None = not probed, False = unavailable


def katex_env():
    """Environment in which node can resolve katex, or False if it cannot."""
    global _katex_env
    if _katex_env is not None:
        return _katex_env
    paths = [VENDOR]
    root = subprocess.run(['npm', 'root', '-g'], capture_output=True, text=True)
    if root.returncode == 0:
        paths.append(root.stdout.strip())
    env = dict(os.environ)
    env['NODE_PATH'] = os.pathsep.join(
        [p for p in paths if p] + ([env['NODE_PATH']] if 'NODE_PATH' in env else []))
    probe = subprocess.run(['node', '-e', 'require.resolve("katex")'],
                           capture_output=True, text=True, env=env)
    if probe.returncode != 0:
        print('note: skipping --katex; install the module with\n'
              f'      mkdir -p {os.path.dirname(VENDOR)} && '
              f'cd {os.path.dirname(VENDOR)} && npm install katex',
              file=sys.stderr)
        _katex_env = False
    else:
        _katex_env = env
    return _katex_env


def katex_check(path, exprs):
    env = katex_env()
    if env is False:
        return []
    payload = json.dumps([{'line': l, 'tex': unescape(t), 'display': d}
                          for l, t, d in exprs])
    res = subprocess.run(['node', '-e', NODE_SCRIPT, payload],
                         capture_output=True, text=True, env=env)
    if res.returncode != 0:
        print(f'note: --katex failed: {res.stderr.strip()[:200]}', file=sys.stderr)
        return []
    return [Finding(path, e['line'], 'katex', e['msg'], '', 'fix the LaTeX')
            for e in json.loads(res.stdout)]


def api_check(path, raw):
    """Ground truth: does GitHub's Markdown split any math span with emphasis?"""
    masked = re.sub(r'```.*?```', 'CODEBLOCK', raw, flags=re.S)
    masked = re.sub(r'^\$\$\n.*?\n\$\$', 'DISPLAYMATH', masked, flags=re.S | re.M)
    req = urllib.request.Request(
        'https://api.github.com/markdown/raw', data=masked.encode('utf-8'),
        headers={'Content-Type': 'text/plain', 'User-Agent': 'mdmath-check'})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            out = r.read().decode('utf-8')
    except Exception as exc:
        print(f'note: skipping --api ({exc})', file=sys.stderr)
        return []
    found = []
    for m in re.finditer(r'<(em|strong)>(.*?)</\1>', out, re.S):
        body = m.group(2)
        if any(c in body for c in '$\\{}^'):
            found.append(Finding(
                path, 0, 'emphasis/api',
                'GitHub split a math span with emphasis',
                html.unescape(body)[:100].replace('\n', ' '),
                'keep every "_" in inline math after a letter or digit'))
    return found


def main():
    ap = argparse.ArgumentParser(
        description='Find LaTeX math that GitHub Markdown renders incorrectly.',
        epilog='See the header of this script for what each check means.')
    ap.add_argument('files', nargs='+', metavar='FILE.md')
    ap.add_argument('--katex', action='store_true',
                    help='also parse each expression with KaTeX (needs node + katex)')
    ap.add_argument('--api', action='store_true',
                    help='also verify against GitHub\'s Markdown API (needs network)')
    ap.add_argument('-q', '--quiet', action='store_true', help='only report problems')
    args = ap.parse_args()

    colour = sys.stdout.isatty()
    total = 0
    for path in args.files:
        try:
            findings, n = check_file(path, args.katex, args.api)
        except OSError as exc:
            print(f'{path}: {exc}', file=sys.stderr)
            total += 1
            continue
        for f in findings:
            f.report(colour)
        total += len(findings)
        if not findings and not args.quiet:
            print(f'ok  {path}  ({n} expressions)')
    if total:
        print(f'\n{total} problem(s) found.')
    return 1 if total else 0


if __name__ == '__main__':
    sys.exit(main())
