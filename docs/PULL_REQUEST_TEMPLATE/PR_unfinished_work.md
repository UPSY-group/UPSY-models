# Unfinished work

<!--
Use this template for PRs that are a step towards a bigger feature, but don't finish it yet.
That's perfectly fine (and encouraged!), as long as it doesn't break anything that already exists.
Text between these arrows won't show up in the PR; feel free to leave it in.
-->

## Which project is this part of?

<!-- Link to the GitHub Discussion where you described the feature you're building. -->

## What does this PR add?

<!-- What does this step do, and how does it fit into the bigger picture?
     E.g. "Adds the routine that calculates polygon intersections; it isn't called anywhere yet,
     but will be used for the conservative remapping in the next PR." -->

## What's next?

<!-- What is left to do after this PR? A short list is enough. -->

## Is the new code switched off?

<!-- Tick one: -->

- [ ] The new code isn't called anywhere yet.
- [ ] The new code is called, but switched off by default through a config parameter: `C%...`

## Checklist

- [ ] The code compiles, and the test suite passes.
- [ ] Model output is unchanged with the default config: all checksums and test references are the same.
- [ ] This PR contains one step of the project, not several.
- [ ] New code follows the [style guide](https://github.com/UPSY-group/UPSY-models/blob/main/docs/STYLE_GUIDE.md).
- [ ] Placeholders that crash on purpose (e.g. `call crash('fixme')`) cannot be reached with the default config.
