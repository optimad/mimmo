<h1><a class="anchor" id="git"></a>
Git Repository Practices</h1>
<p>As most of our code repositories uses git as the revision control system, it is important to decide on a workflow that can be followed by the individual developer. The way that any individual developer interact with the upstream git repository can have an important impact on other developers and the ability to identify and manage individual changes. This set of guidelines and practices attempts to establish some standards for how developers will interact with the upstream git repository.</p>
<h2><a class="anchor" id="commits"></a>
Making Repository Commits</h2>
<p>As a general rule, developers should update frequently, and commit changes often. However, the repository should always remain in a state where the code can be compiled. Most of the time, the code should also successfully execute "make check" run from the top-level directory. If you commit code that violates this principal, it should be your first priority to return the repository code to a compilable state, and your second priority to make sure "make check" runs without errors. Although it would be possible and many software projects do it, we prefer not to force successful execution of the test suite before every commit. Developers should make every effort to avoid having to impose this constraint, by running a make check before every commit.</p>
<p>Commits to the repository should also come with a non-trivial and useful log message.</p>
<h2><a class="anchor" id="outside-master"></a>
Working Outside the Master Branch</h2>
<p>A critical concept is that all changes shall be developed outside of the master<sup>1</sup> branch. Whether they are in a different branch of the upstream<sup>2</sup> repository (gitflow) or a branch of an entirely different fork (forking workflow) is secondary. This is a well-established concept regardless of the workflow being adopted, and allows a number of other benefits as described below.</p>
<h3><a class="anchor" id="fork"></a>
Working on a Different Fork</h3>
<p>There are a number of benefits of working on a different fork rather than a branch of the upstream repo, although not strictly technical:</p><ul>
<li>Developers, particularly new developers, are liberated from the constant oversight of others as they explore new code options. The impact of this may depend on an individual developer’s personality, but for some it creates a refuge where they can be more free and creative.</li>
<li>Similarly, assuming that all changesets in the upstream repo are communicated to the entire development team, the team is spared a noisy stream of notifications and can focus their attention on the rarer occurrence of a pull request notification.</li>
</ul>
<h3><a class="anchor" id="pr"></a>
All Changes are Committed by Pull Request</h3>
<p>Although this can be imposed technically by limiting the authority to alter the upstream repo (as in the forking workflow), a healthy developer community can also simply rely on convention. The advantage of doing it by convention rather than by restriction is that it is easier to distribute the load of reviewing and accepting changes. A critical consequence of this decision is that all code is reviewed before it is committed to the upstream master branch. This has benefits to overall quality in two related ways:</p><ul>
<li>the code under review will improve due to the review itself, and</li>
<li>those involved in the review will maintain a broad awareness of the code base resulting in better contributions from them.</li>
</ul>
<p>This practice does, however, place a substantial burden on the developers to perform timely reviews of the pull requested (PR’ed) code. PR’s that languish without sufficient review have a number of negative consequences:</p><ul>
<li>they need to be refreshed simply to keep them up-to-date with the possibly advancing upstream/master</li>
<li>they may delay further development on similar or related features</li>
<li>they breed frustration in the original developer, undermining the community as a whole. github provides powerful collaboration tools that greatly facilitate this process.</li>
</ul>
<p><sup>1</sup> Although a repository may choose a different name for its main development branch, this document will refer to that as the “master” branch.</p>
<p><sup>2</sup> For this discussion, the “upstream” repo will refer to the centralized authoritative repository used to synchronize changes.</p>
<h2><a class="anchor" id="git-mechanics"></a>
Some Git Mechanics to Keep it Clean</h2>
<p>Given the above practices, there are some mechanical details that can help ensure that the upstream/master repository is always in a state that facilitates all repository actions and interactions.</p>
<ol type="1">
<li><p class="startli">Feature branches being used for development should be kept up-to-date with the upstream/master by rebase only. When a feature branch is rebased against the upstream/master, all changes in the upstream/master are inserted into the feature branch at a point in its history that is prior to any of the changes of the feature branch. This can require conflict resultion as the feature branch changes are “replayed” on top of the new upstream/master in its current state. The primary advantage of this policy is that it keeps all of the feature branch changes contiguous. If, by contrast, the upstream/master is merged into the feature branch, the recent changes in the upstream/master become woven into the set of changes in the feature branch. This can make it more difficult to isolate changes later on.</p>
<p class="startli">Strict adoption of this practice is important since a single merge into a feature branch that is then merged back into the upstream/master can make it nearly impossible for others to rebase.</p>
<p class="startli">A typical workflow with pull-request might look like this, all using the command-line, except for submitting the final pull request. Note that there is never a merge operation.</p><ol type="a">
<li>synchronize your local <code>master</code> branch before anything else <div class="fragment"><div class="line">$ git checkout master</div><div class="line">$ git fetch upstream</div><div class="line">$ git rebase upstream/master</div></div><!-- fragment --></li>
<li>now create a new feature branch from master <div class="fragment"><div class="line">$ git checkout -b my_feature_branch master</div></div><!-- fragment --></li>
<li>now make changes, editing A.cpp, B.hpp, C.cpp</li>
<li>now add/commit your changes to your local feature branch <div class="fragment"><div class="line">$ git add A.cpp B.hpp C.cpp</div><div class="line">$ git commit -m “Make sure you have a good commit message”</div></div><!-- fragment --></li>
<li>push your changes to your feature branch on your fork (often called <code>origin</code>) <div class="fragment"><div class="line">$ git push origin my_feature_branch</div></div><!-- fragment --></li>
<li>make more changes, editing B.hpp, D.hpp, E.cpp</li>
<li>add/commit your changes to your local feature branch <div class="fragment"><div class="line">$ git add B.hpp D.hpp E.cpp</div><div class="line">$ git commit -m “Be sure you have another good commit message”</div></div><!-- fragment --></li>
<li>push your changes to your freature branch on your fork (often called <code>origin</code>) <div class="fragment"><div class="line">$ git push origin my_feature_ranch</div></div><!-- fragment --></li>
<li><p class="startli">When you are ready to submit a pull request, be sure that your feature branch is up-to-date. This first step may seem redundant but is here to be clear which branch we are acting on </p><div class="fragment"><div class="line">$ git checkout my_feature_branch</div><div class="line">$ git fetch upstream</div><div class="line">$ git rebase upstream/master</div></div><!-- fragment --><p> This may generate conflicts that can be addressed at this point.</p>
<p class="startli">NOTE: This step can be performed at any time and should be performed as often as practical to reduce the scope of potential conflicts.</p>
</li>
<li>push your updated feature branch on your fork (often called <code>origin</code>) <div class="fragment"><div class="line">$ git push origin my_feature_branch</div></div><!-- fragment --> This may require the ‘-f’ option to force the push. (It is frequently necessary to force this push because the act of rebasing will “replay” the commits from the feature branch on top of the master, leading to different commit hashes. Each of the commits will contain the same actual information, but because it has a different set of hashes, git will think there is an inconsistency and ask you to force the change.)</li>
<li>Submit a pull request on github, from your fork to the fathomteam fork.</li>
</ol>
</li>
<li><p class="startli">When ready to be adopted into the upstream/master, feature branches should be combined by merge only. This adds the changeset to the end of the upstream/master as a set of individual commits but in a contiguous block.</p>
<p class="startli">A typical workflow to merge a pull-request might look like this, all using the command-line.</p><ol type="a">
<li>synchronize your local <code>master</code> branch before anything else (just because it’s never a bad idea!) <div class="fragment"><div class="line">$ git checkout master</div><div class="line">$ git fetch upstream</div><div class="line">$ git rebase upstream/master</div></div><!-- fragment --></li>
<li>add a remote for the user with the pull-request, perhaps the user is ‘other_user’ <div class="fragment"><div class="line">$ git remote add other_user \</div><div class="line">      git@bitbucket.org:other_user/moab.git</div></div><!-- fragment --></li>
<li>fetch the other users repo <div class="fragment"><div class="line">$ git fetch other_user</div></div><!-- fragment --></li>
<li>check out their feature branch <div class="fragment"><div class="line">$ git checkout -b pr_feature_branch \</div><div class="line">    other_user/feature_branch</div></div><!-- fragment --></li>
<li>confirm that it is up-to-date with the master. This first step may seem redundant but is here to be clear which branch we are acting on <div class="fragment"><div class="line">$ git checkout pr_feature_branch</div><div class="line">$ git fetch upstream</div><div class="line">$ git rebase upstream/master</div></div><!-- fragment --> This may generate conflicts that can be addressed at this point. You may want to request the original author (other_user) take care of these.</li>
<li>once confirmed that it’s up-to-date with master, review this branch including: -reading the code -building the code -running tests</li>
<li>once satisfied that the code meets the necessary standards and that all required/requested changes are fully incorporated into other_users’s feature branch, merge it into master <div class="fragment"><div class="line">$ git checkout master</div></div><!-- fragment --> The next two steps may seem redundant but provide some QA <div class="fragment"><div class="line">$ git fetch upstream</div><div class="line">$ git rebase upstream/master</div><div class="line">$ git merge other_user/feature_branch</div></div><!-- fragment --></li>
<li>push those changes into the master branch on bitbucket <div class="fragment"><div class="line">$ git push upstream/master</div></div><!-- fragment --></li>
</ol>
</li>
<li>When a pull request is open for review, any changes to the feature branch will automatically update the pull request. This is the appropriate way for a developer to respond to requests for changes that occur through the PR process.</li>
<li>If a developer has ongoing work that is based on a feature branch that is under consideration in an open PR, a new feature branch (B) should be created that is based on the previous feature branch (A). Moreover, as changes are made to the original feature branch (A) due to the review process, the new feature branch (B) should be kept up-to-date by rebase against feature branch (A). This keeps all subsequent changes of (B) downstream from the changes in (A). Once feature branch (A) has been adopted into the upstream/master, the new feature branch (B) can start being rebased against the upstream/master instead.</li>
<li>When a repo is forked, its branches are not automatically synchronized with the corresponding branches on the upstream repo. This requires a manual process of synchronization via a local clone. Assuming that the local repo’s branch has the same name as the upstream branch (&lt;branch&gt;), and that the fork is known as “origin”: <div class="fragment"><div class="line">$ git fetch upstream</div><div class="line">$ git checkout &lt;branch&gt;</div><div class="line">$ git rebase upstream/&lt;branch&gt;</div><div class="line">$ git push origin &lt;branch&gt;</div></div><!-- fragment --> The decision of which branches to keep up-to-date is up to the developers. Developers may choose to delete some branches from their own fork to avoid (a) the need to update it and (b) accidentally assuming that it is up-to-date.</li>
<li>When rebasing, it is not uncommon to encounter conflicts. This will interrupt the rebasing process, and each conflicted file will contain conflict blocks. You will be offered three choices:<ul>
<li>manually resolve the conflicts, add the resolved files with git add, and then git rebase &ndash;continue (do not commit the resolved files!)</li>
<li>abort the rebasing process entirely with git rebase &ndash;abort</li>
<li>skip the commit that causes the conflict, assuming that you are sure that this is the right thing to do, with git rebase &ndash;skip</li>
</ul>
</li>
<li>github offers a number of buttons/tools to manage changes between branches and forks. Some of these operate in ways that are contrary to the practices recommended here, and others are consistent with these practices. In general, it is best to know how to do all of those operations with the command-line instead of relying on github, as it gives you full control over important details.</li>
<li>During code development, it might be necessary to work on the same branch on different machines. The workflow to update the local branch is to first fetch the remote changes and then perform a hard reset. <div class="fragment"><div class="line">$ git fetch origin</div><div class="line">$ git reset --hard origin/branch_name</div></div><!-- fragment --> One should be careful with the branch name as a hard reset would overwrite all changes in the working directory.</li>
</ol>
